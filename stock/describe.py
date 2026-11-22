"""Shape descriptors of a scrap.

These are the numbers the rest of the pipeline uses to reason about a scrap without testing its
full outline:

- area, perimeter, vertex count;
- inradius and incenter: the largest circle inside the scrap. Any panel whose circumradius is at
  most the inradius fits in the scrap, which is the coarse check behind the size field;
- circumradius: the smallest circle that encloses the scrap;
- length, width, elongation and axis angle: the minimum-area enclosing rectangle, its sides
  (length >= width), their ratio, and the direction of the long side;
- convexity: area divided by convex hull area (1 for a convex scrap).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import shapely
from shapely.geometry import Polygon
from shapely.geometry.polygon import orient

# Stopping tolerance of the inscribed-circle search, relative to sqrt(area). The radius it returns
# is within this distance of the true inradius. Tighter is much slower on scraps with many
# axis-aligned edges, where many centres tie (remnants: 2.7 ms at 1e-4, 25 ms at 1e-6).
INRADIUS_REL_TOL = 1e-4


def to_polygon(coords) -> Polygon:
    """Check a scrap outline and return it as a counter-clockwise shapely Polygon.

    coords: (n, 2) vertices in order, optionally closed (last vertex equal to the first), or a
    shapely Polygon. Repeated consecutive vertices are dropped. Raises ValueError for an empty
    input, holes, non-finite coordinates, fewer than 3 distinct vertices, or an outline that is not
    a simple polygon (self-intersecting, or degenerate with zero area).
    """
    if isinstance(coords, Polygon):
        if coords.is_empty:
            raise ValueError("scrap is empty")
        if len(coords.interiors) > 0:
            raise ValueError("scraps cannot have holes")
        coords = coords.exterior.coords
    pts = np.asarray(coords, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError(f"expected (n, 2) vertices, got shape {pts.shape}")
    if not np.isfinite(pts).all():
        raise ValueError("scrap vertices must be finite")
    # drop each vertex equal to the next one (cyclically); this also removes a closing vertex
    pts = pts[np.any(pts != np.roll(pts, -1, axis=0), axis=1)]
    if len(pts) < 3:
        raise ValueError(f"a scrap needs at least 3 distinct vertices, got {len(pts)}")
    poly = Polygon(pts)
    if not poly.is_valid:
        raise ValueError(f"invalid scrap outline: {shapely.is_valid_reason(poly)}")
    return orient(poly, sign=1.0)


def inscribed_circle(poly: Polygon, tol: float | None = None) -> tuple[np.ndarray, float]:
    """Center and radius of the largest circle inside poly (GEOS maximum inscribed circle).

    tol: search tolerance in model units; defaults to INRADIUS_REL_TOL * sqrt(area). The center of
    a scrap with a flat stretch (e.g. a long rectangle) is not unique; any valid center is returned.
    """
    if tol is None:
        tol = INRADIUS_REL_TOL * np.sqrt(poly.area)
    line = shapely.maximum_inscribed_circle(poly, tolerance=tol)
    return np.asarray(line.coords[0], dtype=float), float(line.length)


def min_area_rectangle(poly: Polygon) -> tuple[float, float, float]:
    """Minimum-area enclosing rectangle of poly: (length, width, angle), with length >= width.

    angle is the direction of the long side, in radians in [0, pi). One side of the optimal
    rectangle lies along an edge of the convex hull (rotating calipers), so trying every hull edge
    direction gives the exact answer. The optimum need not be unique (for a triangle, all three
    edges give the same area), so among rectangles within a relative 1e-9 of the minimum area the
    most elongated one wins.
    """
    hull = np.asarray(poly.convex_hull.exterior.coords, dtype=float)[:-1]
    edges = np.roll(hull, -1, axis=0) - hull
    # a rectangle is unchanged by a 90 degree turn, so directions modulo pi/2 suffice
    angles = np.unique(np.arctan2(edges[:, 1], edges[:, 0]) % (np.pi / 2))
    c, s = np.cos(angles)[:, None], np.sin(angles)[:, None]
    x, y = hull[:, 0][None, :], hull[:, 1][None, :]
    along = x * c + y * s            # hull coordinates along each candidate direction
    across = -x * s + y * c          # and perpendicular to it
    w_along = along.max(axis=1) - along.min(axis=1)
    w_across = across.max(axis=1) - across.min(axis=1)
    area = w_along * w_across
    elongation = np.maximum(w_along, w_across) / np.minimum(w_along, w_across)
    k = int(np.argmax(np.where(area <= area.min() * (1 + 1e-9), elongation, -np.inf)))
    if w_along[k] >= w_across[k]:
        return float(w_along[k]), float(w_across[k]), float(angles[k] % np.pi)
    return float(w_across[k]), float(w_along[k]), float((angles[k] + np.pi / 2) % np.pi)


@dataclass(frozen=True)
class Descriptors:
    area: float
    perimeter: float
    n_vertices: int
    inradius: float
    incenter: tuple[float, float]
    circumradius: float
    length: float
    width: float
    elongation: float
    axis_angle: float
    convexity: float

    def to_dict(self) -> dict:
        """JSON-ready dict with camelCase keys, like the metrics in web/data/*.meta.json."""
        return {
            "area": self.area,
            "perimeter": self.perimeter,
            "nVertices": self.n_vertices,
            "inradius": self.inradius,
            "incenter": list(self.incenter),
            "circumradius": self.circumradius,
            "length": self.length,
            "width": self.width,
            "elongation": self.elongation,
            "axisAngle": self.axis_angle,
            "convexity": self.convexity,
        }


def describe(scrap, inradius_tol: float | None = None) -> Descriptors:
    """Descriptors of one scrap, given as (n, 2) vertices or a shapely Polygon."""
    poly = to_polygon(scrap)
    center, r_in = inscribed_circle(poly, inradius_tol)
    length, width, angle = min_area_rectangle(poly)
    return Descriptors(
        area=float(poly.area),
        perimeter=float(poly.length),
        n_vertices=len(poly.exterior.coords) - 1,
        inradius=r_in,
        incenter=(float(center[0]), float(center[1])),
        circumradius=float(shapely.minimum_bounding_radius(poly)),
        length=length,
        width=width,
        elongation=length / width,
        axis_angle=angle,
        convexity=float(poly.area / poly.convex_hull.area),
    )
