"""Per-face eta for PH-CPF runs: matlab/run_phcpf.m reads a text file with one value per face.

    python -m sizing.eta metro --eta 0.008 --max-aspect 2 --out results/metro_eta_4x_a2.txt
    python -m sizing.eta metro --eta 0.008 --smooth 200 --out results/metro_eta_4x_smooth.txt

On a face the footprint's aspect is ((k_big^2 + eps^2) / (k_small^2 + eps^2))^(1/4), with
eps = delta / eta: it falls as eps grows. A larger eta therefore stretches footprints, and the
smallest eps that keeps the aspect at most A is

    eps_A^2 = (k_big^2 - A^4 k_small^2) / (A^4 - 1)      (zero if the face is already within A).

`limit_aspect` lowers eta where needed so that every footprint stays within A; elsewhere it keeps
eta. Values are fractions of the surface area, like run_phcpf's eta.

From a target panel area to eta (docs/plan.md, Phase 1, item 2). The panel on a face has area
c delta / ((k1^2 + eps^2)(k2^2 + eps^2))^(1/4), with c = HEX_AREA * PANEL_K^2, so the eps that gives
area A solves a quadratic in x = eps^2:

    (k1^2 + x)(k2^2 + x) = (c delta / A)^4.

`eta_for_area` solves it per face. A target larger than the curvature allows has no solution
(x <= 0); there eta stops at ETA_MAX and the panel is as large as the curvature lets it be.

`smooth_eta` smooths a size field before converting it: it takes the panels a given eta makes,
smooths the log of their area with Laplacian steps over the face adjacency, and converts the result
back. The metric's size then follows the smoothed target, not the per-face curvature estimate,
which jumps between neighbouring faces. `log_jump` measures that roughness.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from sizing.footprint import (DEFAULT_DELTA, DEFAULT_ETA, HEX_AREA, PANEL_K, Surface, footprint,
                              load_surface, panel_count)

ETA_MAX = 16 * DEFAULT_ETA     # the largest eta a conversion gives (the 1x-16x ramp's top rounded)


def limit_aspect(surface: Surface, eta=DEFAULT_ETA, max_aspect: float = 2.0,
                 delta: float = DEFAULT_DELTA) -> np.ndarray:
    """Per-face eta, at most `eta` (one value or one per face), lowered where a footprint would be
    more stretched than `max_aspect` (long / short semi-axis) until it is exactly that."""
    if max_aspect <= 1:
        raise ValueError("max_aspect must be greater than 1")
    footprint(surface, delta, eta)                       # validates delta and eta
    eta = np.broadcast_to(np.asarray(eta, dtype=float), (surface.n_faces,))
    d = delta * surface.span
    ks = np.minimum(surface.kmin, surface.kmax)
    kb = np.maximum(surface.kmin, surface.kmax)
    a4 = max_aspect**4
    eps_a = np.sqrt(np.maximum(kb**2 - a4 * ks**2, 0.0) / (a4 - 1))
    with np.errstate(divide="ignore"):
        eta_a = np.where(eps_a > 0, d / (eps_a * surface.area), np.inf)
    return np.minimum(eta, eta_a)


def face_adjacency(faces: np.ndarray) -> np.ndarray:
    """(m, 2) pairs of faces that share an edge."""
    faces = np.asarray(faces)
    e = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    fid = np.tile(np.arange(len(faces)), 3)
    key = e[:, 0] * (int(faces.max()) + 1) + e[:, 1]
    o = np.argsort(key, kind="stable")
    same = key[o][1:] == key[o][:-1]
    return np.sort(np.c_[fid[o][:-1][same], fid[o][1:][same]], axis=1)


def log_jump(values: np.ndarray, adj: np.ndarray) -> np.ndarray:
    """|log v_f - log v_g| across each pair of neighbouring faces."""
    lv = np.log(np.asarray(values, dtype=float))
    return np.abs(lv[adj[:, 0]] - lv[adj[:, 1]])


def smooth_log(values: np.ndarray, adj: np.ndarray, steps: int, lam: float = 0.5) -> np.ndarray:
    """Positive per-face values after `steps` Laplacian steps on their log: each step moves a face's
    log a fraction lam toward the mean over its neighbours."""
    if steps < 0 or not 0 < lam <= 1:
        raise ValueError("steps must be >= 0 and lam in (0, 1]")
    lv = np.log(np.asarray(values, dtype=float))
    n = len(lv)
    deg = np.maximum(np.bincount(adj.ravel(), minlength=n), 1)
    for _ in range(steps):
        total = np.bincount(adj[:, 0], lv[adj[:, 1]], n) + np.bincount(adj[:, 1], lv[adj[:, 0]], n)
        lv = (1 - lam) * lv + lam * total / deg
    return np.exp(lv)


def eta_for_area(surface: Surface, target, delta: float = DEFAULT_DELTA, max_aspect: float | None = None,
                 eta_max: float = ETA_MAX) -> np.ndarray:
    """Per-face eta whose panel has the target area (m^2, one value or one per face).

    Where the target is larger than the curvature allows, or would need eta above eta_max, eta is
    eta_max and the panel falls short of the target. With max_aspect, eta is lowered further where
    the panel would be more stretched than that (as in limit_aspect), which also makes it smaller.
    """
    target = np.broadcast_to(np.asarray(target, dtype=float), (surface.n_faces,))
    if not np.all(target > 0) or eta_max <= 0:
        raise ValueError("target and eta_max must be positive")
    d = delta * surface.span
    ks = np.minimum(surface.kmin, surface.kmax)
    kb = np.maximum(surface.kmin, surface.kmax)
    r = (HEX_AREA * PANEL_K**2 * d / target)**4
    q = np.sqrt((kb**2 - ks**2)**2 + 4 * r)
    x = 2 * (r - kb**2 * ks**2) / (kb**2 + ks**2 + q)        # positive root, without cancellation
    if max_aspect is not None:
        if max_aspect <= 1:
            raise ValueError("max_aspect must be greater than 1")
        a4 = max_aspect**4
        x = np.maximum(x, np.maximum(kb**2 - a4 * ks**2, 0.0) / (a4 - 1))
    x = np.maximum(x, (d / (eta_max * surface.area))**2)
    return d / (np.sqrt(x) * surface.area)


def smooth_eta(surface: Surface, eta=DEFAULT_ETA, steps: int = 200, max_aspect: float | None = None,
               delta: float = DEFAULT_DELTA, eta_max: float = ETA_MAX) -> np.ndarray:
    """Per-face eta for a smoothed version of the panels `eta` makes: the log of their area gets
    `steps` Laplacian steps, and eta_for_area converts it back."""
    target = footprint(surface, delta, eta).panel_area()
    target = smooth_log(target, face_adjacency(surface.faces), steps)
    return eta_for_area(surface, target, delta, max_aspect, eta_max)


def write_eta(path: Path, eta: np.ndarray) -> Path:
    """One value per line, as run_phcpf's 'eta' option reads it."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, np.asarray(eta, dtype=float).ravel(), fmt="%.10g")
    return path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m sizing.eta",
                                 description="Write a per-face eta file for matlab/run_phcpf.m.")
    ap.add_argument("model")
    ap.add_argument("--eta", type=float, default=DEFAULT_ETA, help=f"fraction of the surface area (default {DEFAULT_ETA:g})")
    ap.add_argument("--max-aspect", type=float, help="largest footprint aspect (long / short semi-axis)")
    ap.add_argument("--smooth", type=int, default=0, metavar="STEPS",
                    help="smooth the panel sizes with this many Laplacian steps, then convert back to eta")
    ap.add_argument("--delta", type=float, default=DEFAULT_DELTA, help=f"fraction of the diagonal (default {DEFAULT_DELTA:g})")
    ap.add_argument("--out", type=Path, required=True, help="text file to write")
    args = ap.parse_args(argv)

    surf = load_surface(args.model)
    if args.smooth > 0:
        eta = smooth_eta(surf, args.eta, args.smooth, args.max_aspect, args.delta)
    elif args.max_aspect is not None:
        eta = limit_aspect(surf, args.eta, args.max_aspect, args.delta)
    else:
        eta = np.full(surf.n_faces, args.eta)
    write_eta(args.out, eta)
    fp = footprint(surf, args.delta, eta)
    lowered = eta < args.eta
    jump = log_jump(np.sqrt(fp.a1 * fp.a2), face_adjacency(surf.faces))
    print(f"{args.model}: eta {eta.min():.3g}..{eta.max():.3g} (median {np.median(eta):.3g}), below "
          f"{args.eta:g} on {lowered.mean():.0%} of faces ({surf.face_area[lowered].sum() / surf.area:.0%} "
          f"of the area); about {panel_count(surf, fp):.0f} panels; aspect max {fp.aspect.max():.2f}; "
          f"log-size jump across edges median {np.median(jump):.4f}; {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
