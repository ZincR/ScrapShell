"""Footprints: the panel size PH-CPF aims for on each face, and the panel it then makes.

    python -m sizing.footprint                  # every model at the 15 m span, baseline parameters
    python -m sizing.footprint --eta 0.008      # the same with a larger maximum element area

PH-CPF sizes panels with a per-face metric g, diagonal in the principal frame with entries
kt_i / delta, where kt_i = sqrt(k_i^2 + eps^2) and eps = delta / eta (SHAPEOP.curvature_sizing,
CPF.set_param). The field asks U'gU = 1, so the footprint on a face is the ellipse with
semi-axes a_i = sqrt(delta / kt_i) along the principal directions. delta is a fraction of the
bounding-box diagonal and eta a fraction of the surface area, as in run_phcpf.m (0.5% and 0.2%).

eps caps both semi-axes at sqrt(eta). Where a principal curvature is smaller than eps the cap, not
the curvature, sets the footprint: with the baseline's eta that is most of the surface on most
models (`capped_fraction`). Raising eta lets footprints grow, but it also stretches them: the
aspect sqrt(kt_max / kt_min) rises from 1 toward sqrt(k_max / k_min) as eps falls.

The extension sets eta per face, from the stock (docs/plan.md, Phase 1), so `footprint` accepts one
eta per face as well as a single value. On each face eta then sets both the panel's size and its
stretch, and never makes a panel larger than the curvature allows (a_i <= sqrt(delta / |k_i|)).

The panels the pipeline makes are, to a good approximation, regular hexagons in the metric g with
circumradius PANEL_K: in 3D, a hexagon stretched by a_1 and a_2 along the principal directions.
PANEL_K was fitted on the baseline's planarized panels (`calibrate`): 0.68-0.71 on each of the six
models, and the predicted panel count is within 6% of the baseline's.

Everything is in meters, with the model at a common span (stock/place.py). The curvatures are
recovered from the sizing ellipses the viewer data stores (stages.guide.axis1, axis2), so no
MATLAB run is needed.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import shapely
from scipy.spatial import cKDTree

from stock.model_ref import WEB_DATA, decode_array, flatten_panel
from stock.place import DEFAULT_SPAN, model_ids

DEFAULT_DELTA = 0.005    # PH-CPF delta: fraction of the bounding-box diagonal (runme.m)
DEFAULT_ETA = 0.002      # PH-CPF eta: maximum element area, fraction of the surface area (runme.m)
PANEL_K = 0.70           # circumradius of a panel in the metric g (fitted by calibrate())
HEX_AREA = 1.5 * math.sqrt(3)    # area of a regular hexagon with circumradius 1


@dataclass(frozen=True)
class Surface:
    """A model's input mesh and per-face principal curvatures, in meters at a given span."""
    model: str
    label: str
    span: float                  # m: bounding-box diagonal
    meters_per_unit: float
    vertices: np.ndarray         # (nv, 3) m
    faces: np.ndarray            # (nf, 3) vertex indices
    face_area: np.ndarray        # (nf,) m^2
    dmin: np.ndarray             # (nf, 3) unit principal direction of kmin
    dmax: np.ndarray             # (nf, 3) unit principal direction of kmax
    kmin: np.ndarray             # (nf,) |principal curvature| along dmin, 1/m
    kmax: np.ndarray             # (nf,) |principal curvature| along dmax, 1/m
    region: np.ndarray           # (nf,) 0 planar, 1 elliptic, 2 parabolic, 3 hyperbolic (viewer stage 2)

    @property
    def area(self) -> float:
        return float(self.face_area.sum())

    @property
    def n_faces(self) -> int:
        return len(self.faces)

    @property
    def centroids(self) -> np.ndarray:
        return self.vertices[self.faces].mean(axis=1)


def _read(model: str, data_dir: Path) -> dict:
    path = Path(data_dir) / f"{model}.json"
    if not path.exists():
        raise ValueError(f"no viewer data for model {model!r} at {path}")
    return json.loads(path.read_text())


def load_surface(model: str, span: float = DEFAULT_SPAN, data_dir: Path = WEB_DATA) -> Surface:
    """The model's surface at `span` meters, with curvatures recovered from its sizing ellipses."""
    if span <= 0:
        raise ValueError("span must be positive")
    data = _read(model, data_dir)
    st = data["stages"]
    if "guide" not in st:
        raise ValueError(f"model {model!r} has no sizing ellipses (viewer stage 3)")
    X = decode_array(st["input"]["V"])
    T = decode_array(st["input"]["F"])
    face_area = 0.5 * np.linalg.norm(np.cross(X[T[:, 1]] - X[T[:, 0]], X[T[:, 2]] - X[T[:, 0]]), axis=1)

    # invert a = sqrt(delta / sqrt(k^2 + eps^2)) with the parameters of the baseline run, in model units
    bbd = float(data["bbd"])
    delta = data["params"]["delta"] * bbd
    eps = delta / (data["params"]["eta"] * face_area.sum())
    axes, dirs, curv = [], [], []
    for key in ("axis1", "axis2"):
        ax = decode_array(st["guide"][key])
        a = np.linalg.norm(ax, axis=1)
        kt = delta / a**2
        axes.append(a)
        dirs.append(ax / a[:, None])
        curv.append(np.sqrt(np.maximum(kt**2 - eps**2, 0.0)))

    s = span / bbd
    region = decode_array(st["curvature"]["region"]).ravel() if "curvature" in st else np.ones(len(T), int)
    return Surface(model=model, label=data.get("label", model), span=span, meters_per_unit=s,
                   vertices=X * s, faces=T, face_area=face_area * s**2, dmin=dirs[0], dmax=dirs[1],
                   kmin=curv[0] / s, kmax=curv[1] / s, region=region)


@dataclass(frozen=True)
class Footprint:
    """Per-face sizing ellipses of a surface (PH-CPF's metric g), in meters."""
    a1: np.ndarray               # semi-axis along dmin
    a2: np.ndarray               # semi-axis along dmax
    delta: float                 # m
    eta: float | np.ndarray      # m^2: maximum element area, one value or one per face; sqrt(eta)
                                 # caps both semi-axes
    capped: np.ndarray           # bool: the cap, not the curvature, sets the larger semi-axis

    @property
    def major(self) -> np.ndarray:
        return np.maximum(self.a1, self.a2)

    @property
    def minor(self) -> np.ndarray:
        return np.minimum(self.a1, self.a2)

    @property
    def aspect(self) -> np.ndarray:
        return self.major / self.minor

    def panel_area(self) -> np.ndarray:
        """Area of the panel PH-CPF makes on each face."""
        return HEX_AREA * PANEL_K**2 * self.a1 * self.a2

    def panel_radius(self) -> np.ndarray:
        """Enclosing-circle radius of that panel when a hexagon vertex lies on the major axis (the
        largest it gets over the hexagon's orientations; exact for a round footprint)."""
        return PANEL_K * self.major


def footprint(surface: Surface, delta: float = DEFAULT_DELTA, eta=DEFAULT_ETA) -> Footprint:
    """Sizing ellipses for PH-CPF parameters delta (fraction of the diagonal) and eta (fraction of
    the surface area: one value, or an array with one value per face)."""
    eta = np.asarray(eta, dtype=float)
    if eta.ndim > 1 or (eta.ndim == 1 and len(eta) != surface.n_faces):
        raise ValueError(f"eta must be one value or one per face ({surface.n_faces}), got shape {eta.shape}")
    if delta <= 0 or not np.all(eta > 0):
        raise ValueError("delta and eta must be positive")
    d = delta * surface.span
    e = eta * surface.area
    eps = d / e
    a1 = np.sqrt(d / np.sqrt(surface.kmin**2 + eps**2))
    a2 = np.sqrt(d / np.sqrt(surface.kmax**2 + eps**2))
    capped = np.minimum(surface.kmin, surface.kmax) < eps
    return Footprint(a1=a1, a2=a2, delta=d, eta=float(e) if e.ndim == 0 else e, capped=capped)


def panel_count(surface: Surface, fp: Footprint) -> float:
    """Number of panels PH-CPF would make on the whole surface (before stripping the boundary)."""
    return float(np.sum(surface.face_area / fp.panel_area()))


def capped_fraction(surface: Surface, fp: Footprint) -> float:
    """Share of the surface area where the cap sqrt(eta), not the curvature, sets the footprint."""
    return float(surface.face_area[fp.capped].sum() / surface.area)


# --- calibration against the baseline's panels -------------------------------------------------

@dataclass(frozen=True)
class Calibration:
    model: str
    n_panels: int                # baseline panels
    k: float                     # PANEL_K fitted on this model: sqrt(median(area / (HEX_AREA a1 a2)))
    k_p10: float
    k_p90: float
    covered: float               # share of the surface the panels cover (the dual drops a boundary layer)
    predicted: float             # panel_count at PANEL_K, times `covered`

    def to_dict(self) -> dict:
        return {"model": self.model, "nPanels": self.n_panels, "k": self.k, "kP10": self.k_p10,
                "kP90": self.k_p90, "covered": self.covered, "predicted": self.predicted}


def calibrate(model: str, data_dir: Path = WEB_DATA) -> Calibration | None:
    """Fit PANEL_K on one model's baseline panels (None if the baseline did not planarize it).

    Each panel is compared with the footprint of the face nearest its centroid, at the baseline's
    own delta and eta. Units cancel, so the span does not matter.
    """
    data = _read(model, data_dir)
    st = data["stages"]
    if not st.get("planar", {}).get("frames"):
        return None
    surf = load_surface(model, data_dir=data_dir)
    fp = footprint(surf, data["params"]["delta"], data["params"]["eta"])
    V = decode_array(st["planar"]["frames"][-1]["V"]) * surf.meters_per_unit
    deg = decode_array(st["dual"]["D"]).ravel()
    idx = decode_array(st["dual"]["idx"]).ravel()
    starts = np.concatenate([[0], np.cumsum(deg)])
    panels = [V[idx[a:b]] for a, b in zip(starts[:-1], starts[1:])]
    _, face = cKDTree(surf.centroids).query(np.array([p.mean(axis=0) for p in panels]))
    area = np.array([shapely.Polygon(flatten_panel(p)).area for p in panels])
    k2 = area / (HEX_AREA * fp.a1[face] * fp.a2[face])
    covered = float(area.sum() / surf.area)
    return Calibration(model=model, n_panels=len(panels), k=float(np.sqrt(np.median(k2))),
                       k_p10=float(np.sqrt(np.percentile(k2, 10))),
                       k_p90=float(np.sqrt(np.percentile(k2, 90))), covered=covered,
                       predicted=panel_count(surf, fp) * covered)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m sizing.footprint",
                                 description="Footprints of every model, and PANEL_K against the baseline.")
    ap.add_argument("--span", type=float, default=DEFAULT_SPAN, help=f"bounding-box diagonal in m (default {DEFAULT_SPAN:g})")
    ap.add_argument("--delta", type=float, default=DEFAULT_DELTA, help=f"fraction of the diagonal (default {DEFAULT_DELTA:g})")
    ap.add_argument("--eta", type=float, default=DEFAULT_ETA, help=f"fraction of the surface area (default {DEFAULT_ETA:g})")
    ap.add_argument("--models", nargs="+", help="models (default: all in the manifest)")
    args = ap.parse_args(argv)

    print(f"span {args.span:g} m, delta {args.delta:g}, eta {args.eta:g}")
    print(f"{'model':14s} {'area m^2':>9s} {'cap m':>6s} {'capped':>7s} {'major m p10/50/90':>18s} "
          f"{'aspect p50/90':>14s} {'panel m^2':>10s} {'panels':>7s}")
    for m in args.models or model_ids():
        surf = load_surface(m, args.span)
        fp = footprint(surf, args.delta, args.eta)
        maj = np.percentile(fp.major, [10, 50, 90])
        asp = np.percentile(fp.aspect, [50, 90])
        print(f"{m:14s} {surf.area:9.1f} {math.sqrt(fp.eta):6.3f} {100 * capped_fraction(surf, fp):6.0f}% "
              f"{maj[0]:6.3f}{maj[1]:6.3f}{maj[2]:6.3f} {asp[0]:7.2f}{asp[1]:7.2f} "
              f"{np.median(fp.panel_area()):10.3f} {panel_count(surf, fp):7.0f}")

    print(f"\nPANEL_K = {PANEL_K} against the baseline panels (at the baseline's own delta and eta)")
    print(f"{'model':14s} {'k':>6s} {'k p10-p90':>12s} {'covered':>8s} {'predicted':>10s} {'baseline':>9s}")
    for m in args.models or model_ids():
        c = calibrate(m)
        if c is None:
            print(f"{m:14s} (no baseline panels)")
            continue
        print(f"{m:14s} {c.k:6.3f} {c.k_p10:6.3f}-{c.k_p90:5.3f} {100 * c.covered:7.0f}% "
              f"{c.predicted:10.0f} {c.n_panels:9d}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
