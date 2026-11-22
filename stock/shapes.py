"""Random scrap shapes, in four families.

Each family draws a raw outline. `generate_shape` cleans it, checks it against `Limits`, redraws
until it passes, and returns it in a canonical frame: unit area, incenter at the origin, long axis
along x. Scaling to a target size happens later (stock/sizes.py), so the families only decide the
shape and every family covers the same range of sizes.

- remnant: a sheet with 1-3 rectangles or hexagons cut in from its corners or edges, i.e. what is
  left after cutting parts from a sheet. L, U, T and notched shapes. The default family.
- shard: one cell of a random Voronoi diagram, clipped to a sheet, like a broken tile. Convex.
- blob: a star-shaped polygon with jittered angles and radii. Irregular, often non-convex.
- board: a rectangle with a lognormal aspect ratio, like a plank offcut.
"""
from __future__ import annotations

import math
from collections import Counter
from collections.abc import Callable, Mapping
from dataclasses import dataclass

import numpy as np
import shapely
from shapely import affinity
from shapely.geometry import MultiPoint, Polygon, box

from stock.describe import describe, inscribed_circle, to_polygon


@dataclass(frozen=True)
class Limits:
    """What a usable scrap looks like. Lengths are relative to sqrt(area); angles are in degrees."""
    max_vertices: int = 16
    min_fatness: float = 0.2     # inradius / sqrt(area): 0.5 for a square, 0.22 for a 5:1 board
    min_edge: float = 0.03       # shortest edge
    min_angle: float = 20.0      # every interior angle lies in [min_angle, 360 - min_angle]
    min_feature: float = 0.15    # no part (spike, strip, neck) and no notch (slot) narrower than this


DEFAULT_LIMITS = Limits()
DEFAULT_MIX = {"remnant": 0.6, "shard": 0.15, "blob": 0.15, "board": 0.1}

# Simplification tolerance applied to every raw outline, relative to sqrt(area). It removes
# near-collinear vertices and the slivers that boolean operations leave behind.
SIMPLIFY_TOL = 0.005
_MITRE = {"join_style": "mitre", "mitre_limit": 10.0}   # sharp joins: 1 / sin(10 deg) < 10


# --- checks and canonical frame ---------------------------------------------------------------

def interior_angles(pts: np.ndarray) -> np.ndarray:
    """Interior angle at each vertex of a counter-clockwise polygon, in degrees (reflex > 180)."""
    to_prev = np.roll(pts, 1, axis=0) - pts
    to_next = np.roll(pts, -1, axis=0) - pts
    cross = to_next[:, 0] * to_prev[:, 1] - to_next[:, 1] * to_prev[:, 0]
    dot = (to_next * to_prev).sum(axis=1)
    return np.degrees(np.arctan2(cross, dot)) % 360.0


def check(poly, limits: Limits = DEFAULT_LIMITS) -> str | None:
    """None if poly is a usable scrap, otherwise the reason it is not (a short "kind: detail")."""
    try:
        poly = to_polygon(poly)
    except ValueError as err:
        return f"outline: {err}"
    s = math.sqrt(poly.area)
    pts = np.asarray(poly.exterior.coords)[:-1]
    if len(pts) > limits.max_vertices:
        return f"vertices: {len(pts)} > {limits.max_vertices}"
    edges = np.linalg.norm(np.roll(pts, -1, axis=0) - pts, axis=1)
    if edges.min() < limits.min_edge * s:
        return "edge: too short"
    angles = interior_angles(pts)
    if angles.min() < limits.min_angle or angles.max() > 360.0 - limits.min_angle:
        return "angle: too sharp"
    # Morphological opening removes every part narrower than min_feature, and closing fills every
    # notch narrower than it. With mitre joins both leave wider features exactly as they are, so a
    # usable scrap survives both unchanged.
    w = 0.5 * limits.min_feature * s
    opened = poly.buffer(-w, **_MITRE).buffer(w, **_MITRE)
    if opened.geom_type != "Polygon" or opened.area < 0.99 * poly.area:
        return "feature: too thin"
    closed = poly.buffer(w, **_MITRE).buffer(-w, **_MITRE)
    if closed.area > 1.01 * poly.area:
        return "notch: too narrow"
    if inscribed_circle(poly, tol=1e-3 * s)[1] < limits.min_fatness * s:   # a threshold: coarse is enough
        return "fatness: too thin overall"
    return None


def canonical(poly) -> Polygon:
    """Scale to unit area, move the incenter to the origin and turn the long axis onto x."""
    d = describe(poly)
    p = affinity.translate(to_polygon(poly), -d.incenter[0], -d.incenter[1])
    p = affinity.rotate(p, -d.axis_angle, origin=(0, 0), use_radians=True)
    k = 1.0 / math.sqrt(d.area)
    return to_polygon(affinity.scale(p, k, k, origin=(0, 0)))


def _largest_polygon(geom) -> Polygon | None:
    """The largest polygon in geom (boolean operations can return several pieces), or None."""
    if geom is None or geom.is_empty:
        return None
    parts = [g for g in getattr(geom, "geoms", [geom]) if g.geom_type == "Polygon" and not g.is_empty]
    return max(parts, key=lambda g: g.area) if parts else None


# --- families ---------------------------------------------------------------------------------

def _lognormal_in(rng, median: float, sigma: float, lo: float, hi: float) -> float:
    """Lognormal sample conditioned on [lo, hi] (resampled, so there is no pile-up at the ends)."""
    for _ in range(100):
        x = math.exp(rng.normal(math.log(median), sigma))
        if lo <= x <= hi:
            return x
    return min(max(x, lo), hi)


def board(rng, aspect_median: float = 2.0, aspect_sigma: float = 0.5, aspect_max: float = 5.0) -> Polygon:
    """A rectangle whose aspect ratio (long / short side) is lognormal on [1, aspect_max]."""
    return box(0, 0, _lognormal_in(rng, aspect_median, aspect_sigma, 1.0, aspect_max), 1)


def shard(rng, n_sites: int = 10, aspect_max: float = 1.5) -> Polygon:
    """One random cell of the Voronoi diagram of n_sites uniform points, clipped to the sheet."""
    sheet = box(0, 0, rng.uniform(1, aspect_max), 1)
    x0, y0, x1, y1 = sheet.bounds
    sites = MultiPoint(np.c_[rng.uniform(x0, x1, n_sites), rng.uniform(y0, y1, n_sites)])
    cells = shapely.voronoi_polygons(sites, extend_to=sheet).geoms
    return cells[int(rng.integers(len(cells)))].intersection(sheet)


def blob(rng, n_vertices: tuple[int, int] = (5, 12), irregularity: float = 0.6,
         spikiness: float = 0.25) -> Polygon:
    """A star-shaped polygon: n vertices at jittered, evenly spaced angles with noisy radii.

    irregularity in [0, 1): angular jitter, as a fraction of the even spacing (below 1, vertices
    keep their order, so the outline never crosses itself). spikiness: std of the radius around 1.
    """
    if not 0 <= irregularity < 1:
        raise ValueError("irregularity must be in [0, 1)")
    n = int(rng.integers(n_vertices[0], n_vertices[1] + 1))
    ang = 2 * np.pi * (np.arange(n) + rng.uniform(-irregularity / 2, irregularity / 2, n)) / n
    rad = np.clip(rng.normal(1.0, spikiness, n), 0.3, 1.7)
    return Polygon(np.c_[rad * np.cos(ang), rad * np.sin(ang)])


def _hexagon(center, circumradius: float, rotation: float) -> Polygon:
    a = rotation + np.arange(6) * np.pi / 3
    return Polygon(np.c_[center[0] + circumradius * np.cos(a), center[1] + circumradius * np.sin(a)])


def remnant(rng, aspect_max: float = 2.5, n_cuts: tuple[int, int] = (1, 3), p_corner: float = 0.5,
            p_hex: float = 0.3, corner_cut: tuple[float, float] = (0.3, 0.65),
            notch_width: tuple[float, float] = (0.2, 0.5), notch_depth: tuple[float, float] = (0.2, 0.6),
            hex_size: tuple[float, float] = (0.25, 0.55), min_keep: float = 0.35) -> Polygon | None:
    """A W x 1 sheet (W uniform in [1, aspect_max]) with parts cut in from its boundary.

    Each cut is, with probability p_hex, a regular hexagon centred on a random boundary point
    (circumradius hex_size x the short side, random rotation). Otherwise it is a rectangle: from a
    corner with probability p_corner (corner_cut x each side, which leaves an L), or else a notch
    in one edge (notch_width x that edge, notch_depth x the other side, at a random position).
    Every cut crosses the boundary, so no holes form. The largest leftover piece is returned, or
    None if less than min_keep of the sheet is left.
    """
    W, H = rng.uniform(1, aspect_max), 1.0
    sheet = box(0, 0, W, H)
    piece = sheet
    out = 1.0   # rectangle cuts extend this far past the sheet, so they never leave slivers at its edge
    for _ in range(int(rng.integers(n_cuts[0], n_cuts[1] + 1))):
        if rng.random() < p_hex:
            # a boundary point, uniform along the perimeter
            t = rng.uniform(0, 2 * (W + H))
            if t < W:
                c = (t, 0.0)
            elif t < W + H:
                c = (W, t - W)
            elif t < 2 * W + H:
                c = (2 * W + H - t, H)
            else:
                c = (0.0, 2 * (W + H) - t)
            cut = _hexagon(c, rng.uniform(*hex_size) * H, rng.uniform(0, np.pi / 3))
        elif rng.random() < p_corner:
            w, h = rng.uniform(*corner_cut) * W, rng.uniform(*corner_cut) * H
            cx, cy = int(rng.integers(2)), int(rng.integers(2))      # which corner
            x = (-out, w) if cx == 0 else (W - w, W + out)
            y = (-out, h) if cy == 0 else (H - h, H + out)
            cut = box(x[0], y[0], x[1], y[1])
        else:
            side = int(rng.choice(4, p=np.array([W, H, W, H]) / (2 * (W + H))))   # bottom, right, top, left
            along, across = (W, H) if side in (0, 2) else (H, W)
            w, d = rng.uniform(*notch_width) * along, rng.uniform(*notch_depth) * across
            a0 = rng.uniform(-w / 2, along - w / 2)
            if side == 0:
                cut = box(a0, -out, a0 + w, d)
            elif side == 2:
                cut = box(a0, H - d, a0 + w, H + out)
            elif side == 1:
                cut = box(W - d, a0, W + out, a0 + w)
            else:
                cut = box(-out, a0, d, a0 + w)
        piece = _largest_polygon(piece.difference(cut))
        if piece is None:
            return None
    return piece if piece.area >= min_keep * sheet.area else None


FAMILIES: dict[str, Callable[..., Polygon | None]] = {
    "remnant": remnant, "shard": shard, "blob": blob, "board": board,
}


# --- sampling ---------------------------------------------------------------------------------

def generate_shape(rng: np.random.Generator, family: str = "remnant", params: Mapping | None = None,
                   limits: Limits = DEFAULT_LIMITS, max_tries: int = 200,
                   stats: Counter | None = None) -> Polygon:
    """Draw one usable scrap of the given family, in the canonical frame (unit area).

    params: keyword arguments for the family function. stats, if given, counts the outcome of every
    draw ("ok" or the rejection reason), which shows how restrictive the limits are.
    """
    if family not in FAMILIES:
        raise ValueError(f"unknown family {family!r}; choose from {sorted(FAMILIES)}")
    draw = FAMILIES[family]
    for _ in range(max_tries):
        raw = _largest_polygon(draw(rng, **(params or {})))
        if raw is None:
            reason = "draw: nothing left"
        else:
            raw = raw.simplify(SIMPLIFY_TOL * math.sqrt(raw.area), preserve_topology=True)
            reason = check(raw, limits)
        if stats is not None:
            stats[(reason or "ok").split(":")[0]] += 1
        if reason is None:
            return canonical(raw)
    raise RuntimeError(f"no usable {family} scrap in {max_tries} draws; loosen the limits or the "
                       f"family parameters")


def sample_shapes(rng: np.random.Generator, n: int, mix: Mapping[str, float] = DEFAULT_MIX,
                  params: Mapping[str, Mapping] | None = None, limits: Limits = DEFAULT_LIMITS,
                  stats: Counter | None = None) -> list[tuple[str, Polygon]]:
    """n scraps in the canonical frame, each from a family drawn with the weights in mix."""
    names = list(mix)
    unknown = set(names) - set(FAMILIES)
    if unknown:
        raise ValueError(f"unknown families {sorted(unknown)}; choose from {sorted(FAMILIES)}")
    w = np.array([mix[k] for k in names], dtype=float)
    if (w < 0).any() or w.sum() <= 0:
        raise ValueError("mix weights must be non-negative and not all zero")
    picks = rng.choice(len(names), size=n, p=w / w.sum())
    params = params or {}
    return [(names[i], generate_shape(rng, names[i], params.get(names[i]), limits, stats=stats))
            for i in picks]
