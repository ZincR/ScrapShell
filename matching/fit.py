"""Fit test: does a flat panel fit inside a scrap? (docs/plan.md, Phase 2, item 5.)

A panel P fits scrap S if, for some rotation, a translation t puts P + t inside S. The set of such t
is S shrunk by P,

    S (-) P = {t : P + t inside S} = S minus (boundary of S (+) -P),

because P + t, which is connected and contains t, lies inside S exactly when t is in S and P + t
does not meet the boundary. The Minkowski sum of each boundary edge [a, b] with -P is the convex
hull of (a - P) and (b - P) for a convex P; a non-convex panel (a bow-tie) is split into triangles
first. Translation is exact; rotation is sampled (ROTATIONS angles over a full turn).

Pairs are filtered first:
- cannot fit: the scrap is smaller than the panel in area, inradius, circumradius or minimum width
  (each of these only grows from a shape to a shape containing it);
- sure fit: the scrap's inscribed circle can hold the panel's enclosing circle.
Only the pairs in between get the exact test.

Panels come from the viewer data, flattened onto their best-fit plane (stock/model_ref.py) and
scaled to meters with the model at the deposit's 15 m span (stock/place.py).
"""
from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import Polygon

from stock.model_ref import ROOT, load_model_ref
from stock.place import DEFAULT_SPAN, placement

ROTATIONS = 36
DEPOSIT = ROOT / "data" / "stock" / "deposit_s0.json"


def min_width(points) -> float:
    """Smallest distance between two parallel lines enclosing the points (rotating calipers on the
    hull: the narrowest strip is flush with a hull edge)."""
    hull = shapely.MultiPoint(np.asarray(points, dtype=float)).convex_hull
    if hull.geom_type != "Polygon":
        return 0.0
    h = np.asarray(hull.exterior.coords)[:-1]
    e = np.roll(h, -1, axis=0) - h
    n = np.c_[-e[:, 1], e[:, 0]] / np.linalg.norm(e, axis=1)[:, None]
    dist = np.abs(np.einsum("ijk,ik->ij", h[None, :, :] - h[:, None, :], n))
    return float(dist.max(axis=1).min())


@dataclass(frozen=True)
class Shape:
    """Rotation-invariant descriptors used by the filters."""
    area: float
    inradius: float
    circumradius: float
    width: float


def describe_panel(P: np.ndarray) -> Shape:
    poly = Polygon(P)
    return Shape(area=float(poly.area),
                 inradius=float(shapely.maximum_inscribed_circle(poly).length),
                 circumradius=float(shapely.minimum_bounding_radius(poly)),
                 width=min_width(P))


def centered(P: np.ndarray) -> np.ndarray:
    """The panel with its reference point at the origin: the vertex mean, or a point inside the
    panel if the mean falls outside it."""
    P = np.asarray(P, dtype=float)
    c = P.mean(axis=0)
    poly = Polygon(P)
    if not poly.contains(shapely.Point(c)):
        c = np.asarray(poly.representative_point().coords[0])
    return P - c


def _pieces(P: np.ndarray) -> list[np.ndarray]:
    """Convex pieces of a panel: itself if convex, else its triangles."""
    poly = Polygon(P)
    if poly.convex_hull.area - poly.area <= 1e-9 * poly.area:
        return [P]
    tris = shapely.constrained_delaunay_triangles(poly)
    return [np.asarray(t.exterior.coords)[:-1] for t in tris.geoms]


def inner_fit_region(scrap: Polygon, P: np.ndarray):
    """Translations t that put P + t inside the scrap (P around its reference point)."""
    E = np.asarray(scrap.exterior.coords)
    a, b = E[:-1], E[1:]
    hulls = []
    for Q in _pieces(P):
        pts = np.concatenate([a[:, None, :] - Q[None], b[:, None, :] - Q[None]], axis=1)
        hulls.append(shapely.convex_hull(shapely.multipoints(pts)))
    return scrap.difference(shapely.union_all(np.concatenate(hulls)))


def rotate(P: np.ndarray, angle: float) -> np.ndarray:
    c, s = np.cos(angle), np.sin(angle)
    return P @ np.array([[c, s], [-s, c]])


def fit(scrap: Polygon, P: np.ndarray, rotations: int = ROTATIONS):
    """(angle, t) of a placement of the panel P (around its reference point) inside the scrap,
    trying `rotations` angles in turn, or None if none of them fits."""
    for k in range(rotations):
        angle = 2 * np.pi * k / rotations
        region = inner_fit_region(scrap, rotate(P, angle))
        if region.area > 0:
            t = np.asarray(region.representative_point().coords[0])
            return angle, t
    return None


def load_deposit(path: Path = DEPOSIT) -> dict:
    return json.loads(Path(path).read_text())


def model_panels(model: str, span: float = DEFAULT_SPAN) -> list[np.ndarray]:
    """The model's planarized panels, flattened, in meters with the model at `span`."""
    s = placement(model, span).meters_per_unit
    return [p * s for p in load_model_ref(model).panels]


@dataclass
class FitGraph:
    """Which scraps each panel fits: placements[i] maps scrap index -> (angle, t)."""
    n_panels: int
    n_scraps: int
    placements: list[dict]
    n_sure: int
    n_tested: int
    n_tested_fit: int

    def fits(self) -> np.ndarray:
        """(n_panels, n_scraps) bool."""
        m = np.zeros((self.n_panels, self.n_scraps), bool)
        for i, pl in enumerate(self.placements):
            m[i, list(pl)] = True
        return m


_SCRAPS: list = []


def _init(scraps):
    global _SCRAPS
    _SCRAPS = scraps


def _panel_job(args):
    i, P, candidates, sure, rotations = args
    scraps = _SCRAPS
    out, n_fit = {}, 0
    P0 = centered(P)
    c0 = np.asarray(shapely.minimum_bounding_circle(Polygon(P0)).centroid.coords[0])
    for j in sure:                                   # enclosing circle centred on the incenter
        out[int(j)] = (0.0, np.asarray(scraps[j]["incenter"]) - c0)
    for j in candidates:
        r = fit(Polygon(scraps[j]["polygon"]), P0, rotations)
        if r is not None:
            out[int(j)] = r
            n_fit += 1
    return i, out, n_fit


def fit_graph(panels: list[np.ndarray], deposit: dict, rotations: int = ROTATIONS,
              workers: int | None = None) -> FitGraph:
    """Test every panel against every scrap (filters first, exact test for the rest)."""
    scraps = deposit["scraps"]
    S_area = np.array([s["area"] for s in scraps])
    S_in = np.array([s["inradius"] for s in scraps])
    S_circ = np.array([s["circumradius"] for s in scraps])
    S_w = np.array([min_width(s["polygon"]) for s in scraps])
    shapes = [describe_panel(P) for P in panels]
    jobs, n_sure, n_tested = [], 0, 0
    for i, (P, d) in enumerate(zip(panels, shapes)):
        possible = (S_area >= d.area) & (S_in >= d.inradius) & (S_circ >= d.circumradius) & (S_w >= d.width)
        sure = possible & (S_in >= d.circumradius)
        test = possible & ~sure
        n_sure += int(sure.sum())
        n_tested += int(test.sum())
        jobs.append((i, P, np.flatnonzero(test), np.flatnonzero(sure), rotations))
    placements = [None] * len(panels)
    n_fit = 0
    workers = workers or max(1, (os.cpu_count() or 2) - 2)
    if workers == 1:
        _init(scraps)
        results = map(_panel_job, jobs)
    else:
        pool = ProcessPoolExecutor(workers, initializer=_init, initargs=(scraps,))
        results = pool.map(_panel_job, jobs, chunksize=4)
    for i, out, k in results:
        placements[i] = out
        n_fit += k
    if workers != 1:
        pool.shutdown()
    return FitGraph(n_panels=len(panels), n_scraps=len(scraps), placements=placements,
                    n_sure=n_sure, n_tested=n_tested, n_tested_fit=n_fit)
