"""Tests for stock/describe.py: exact values on known shapes, invariances, and input checks."""
import json
import math
import warnings

import numpy as np
import pytest
import shapely
from shapely.geometry import Point, Polygon

from stock.describe import describe, inscribed_circle, min_area_rectangle, to_polygon

SQRT2, SQRT3, SQRT5 = math.sqrt(2), math.sqrt(3), math.sqrt(5)
TOL = 1e-6

SQUARE = [(0, 0), (1, 0), (1, 1), (0, 1)]
RECT_2x1 = [(0, 0), (2, 0), (2, 1), (0, 1)]
HEXAGON = [(math.cos(k * math.pi / 3), math.sin(k * math.pi / 3)) for k in range(6)]  # circumradius 1
TRIANGLE = [(0, 0), (1, 0), (0.5, SQRT3 / 2)]                                          # side 1
L_SHAPE = [(0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2)]                           # 2x2 minus a 1x1 corner


def star_polygon(seed, n=9):
    """Random non-convex polygon, star-shaped about the origin (angle gaps stay below pi)."""
    rng = np.random.default_rng(seed)
    ang = 2 * np.pi * (np.arange(n) + rng.uniform(-0.4, 0.4, n)) / n
    rad = rng.uniform(0.4, 1.0, n)
    return np.c_[rad * np.cos(ang), rad * np.sin(ang)]


def angle_diff(a, b):
    """Distance between two axis directions, modulo pi."""
    d = (a - b) % math.pi
    return min(d, math.pi - d)


# --- known shapes -----------------------------------------------------------------------------

def test_unit_square():
    d = describe(SQUARE)
    assert d.area == pytest.approx(1)
    assert d.perimeter == pytest.approx(4)
    assert d.n_vertices == 4
    assert d.inradius == pytest.approx(0.5, abs=TOL)
    assert d.incenter == pytest.approx((0.5, 0.5), abs=TOL)
    assert d.circumradius == pytest.approx(SQRT2 / 2)
    assert (d.length, d.width, d.elongation) == pytest.approx((1, 1, 1))
    assert d.convexity == pytest.approx(1)


def test_rectangle_2x1():
    d = describe(RECT_2x1)
    assert d.area == pytest.approx(2)
    assert d.perimeter == pytest.approx(6)
    assert d.inradius == pytest.approx(0.5, abs=TOL)
    # the center is not unique: anywhere on the midline between x = 0.5 and x = 1.5
    assert d.incenter[1] == pytest.approx(0.5, abs=TOL)
    assert 0.5 - TOL <= d.incenter[0] <= 1.5 + TOL
    assert d.circumradius == pytest.approx(SQRT5 / 2)
    assert (d.length, d.width, d.elongation) == pytest.approx((2, 1, 2))
    assert angle_diff(d.axis_angle, 0) == pytest.approx(0, abs=1e-12)
    assert d.convexity == pytest.approx(1)


def test_regular_hexagon():
    d = describe(HEXAGON)
    assert d.area == pytest.approx(3 * SQRT3 / 2)
    assert d.perimeter == pytest.approx(6)
    assert d.n_vertices == 6
    assert d.inradius == pytest.approx(SQRT3 / 2, abs=TOL)
    assert d.incenter == pytest.approx((0, 0), abs=TOL)
    assert d.circumradius == pytest.approx(1)
    assert (d.length, d.width) == pytest.approx((2, SQRT3))
    assert d.elongation == pytest.approx(2 / SQRT3)
    assert d.convexity == pytest.approx(1)


def test_equilateral_triangle():
    d = describe(TRIANGLE)
    assert d.area == pytest.approx(SQRT3 / 4)
    assert d.inradius == pytest.approx(SQRT3 / 6, abs=TOL)
    assert d.incenter == pytest.approx((0.5, SQRT3 / 6), abs=TOL)
    assert d.circumradius == pytest.approx(1 / SQRT3)
    assert (d.length, d.width) == pytest.approx((1, SQRT3 / 2))


def test_l_shape():
    d = describe(L_SHAPE)
    assert d.area == pytest.approx(3)
    assert d.perimeter == pytest.approx(8)
    assert d.n_vertices == 6
    # the largest circle sits in the corner square and touches the reflex vertex (1, 1):
    # c = sqrt(2) (1 - c)  ->  c = 2 - sqrt(2)
    assert d.inradius == pytest.approx(2 - SQRT2, abs=TOL)
    assert d.incenter == pytest.approx((2 - SQRT2, 2 - SQRT2), abs=TOL)
    assert d.circumradius == pytest.approx(SQRT2)         # circle through (2, 0) and (0, 2)
    assert (d.length, d.width) == pytest.approx((2, 2))
    assert d.convexity == pytest.approx(3 / 3.5)          # hull cuts the corner triangle of area 1/2


# --- invariances and cross-checks -------------------------------------------------------------

@pytest.mark.parametrize("seed", range(5))
def test_similarity_invariance(seed):
    """Moving, turning and scaling a scrap scales lengths by k and areas by k^2."""
    pts = star_polygon(seed)
    theta, k, t = 0.7 + seed, 3.5, np.array([10.0, -4.0])
    R = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
    a, b = describe(pts), describe(k * pts @ R.T + t)
    assert b.area == pytest.approx(k**2 * a.area)
    for f in ("perimeter", "inradius", "circumradius", "length", "width"):
        assert getattr(b, f) == pytest.approx(k * getattr(a, f), rel=1e-5), f
    assert b.elongation == pytest.approx(a.elongation, rel=1e-5)
    assert b.convexity == pytest.approx(a.convexity)
    assert angle_diff(b.axis_angle, a.axis_angle + theta) == pytest.approx(0, abs=1e-6)
    assert b.incenter == pytest.approx(tuple(k * R @ np.array(a.incenter) + t), abs=1e-4 * k)


@pytest.mark.parametrize("seed", range(5))
def test_inradius_against_brute_force(seed):
    """No grid point inside the scrap is farther from the boundary than the reported inradius,
    and the reported circle lies inside the scrap."""
    poly = to_polygon(star_polygon(seed))
    center, r = inscribed_circle(poly)
    x0, y0, x1, y1 = poly.bounds
    h = (x1 - x0) / 300
    gx, gy = np.meshgrid(np.arange(x0, x1, h), np.arange(y0, y1, h))
    pts = shapely.points(gx.ravel(), gy.ravel())
    inside = shapely.contains(poly, pts)
    best = shapely.distance(pts[inside], poly.exterior).max()
    assert best <= r + 1e-9
    assert r - best <= h                                   # the grid gets within one cell of it
    assert poly.exterior.distance(Point(center)) == pytest.approx(r, abs=1e-6)
    assert poly.contains(Point(center))


@pytest.mark.parametrize("seed", range(5))
def test_min_area_rectangle_matches_shapely(seed):
    poly = to_polygon(star_polygon(seed))
    length, width, _ = min_area_rectangle(poly)
    with warnings.catch_warnings(), np.errstate(all="ignore"):
        warnings.simplefilter("ignore", RuntimeWarning)
        ref = shapely.oriented_envelope(poly).area
    assert length * width == pytest.approx(ref, rel=1e-9)
    assert length >= width


# --- input handling ---------------------------------------------------------------------------

def test_clockwise_input_is_reoriented():
    poly = to_polygon(SQUARE[::-1])
    assert shapely.is_ccw(poly.exterior)
    assert describe(SQUARE[::-1]).area == pytest.approx(1)


def test_closed_ring_and_polygon_inputs_agree():
    a = describe(L_SHAPE)
    assert describe(L_SHAPE + [L_SHAPE[0]]) == a
    assert describe(Polygon(L_SHAPE)) == a


def test_repeated_vertices_are_dropped():
    d = describe([(0, 0), (1, 0), (1, 0), (0, 1), (0, 0)])
    assert d.n_vertices == 3
    assert d.area == pytest.approx(0.5)


@pytest.mark.parametrize("bad, message", [
    ([(0, 0), (1, 0)], "at least 3"),
    ([(0, 0), (1, 0), (1, 0), (0, 0)], "at least 3"),           # only 2 distinct vertices
    ([(0, 0), (1, 0), (2, 0)], "invalid"),                      # collinear: zero area
    ([(0, 0), (1, 1), (1, 0), (0, 1)], "invalid"),              # symmetric bow-tie (signed area 0)
    ([(0, 0), (2, 2), (2, 0), (0, 1)], "invalid"),              # asymmetric bow-tie
    ([(0, 0), (1, 0), (math.nan, 1)], "finite"),
    ([(0, 0, 0), (1, 0, 0), (0, 1, 0)], "shape"),
])
def test_invalid_outlines_are_rejected(bad, message):
    with pytest.raises(ValueError, match=message):
        to_polygon(bad)


def test_holes_are_rejected():
    ring = Polygon([(0, 0), (4, 0), (4, 4), (0, 4)], holes=[[(1, 1), (2, 1), (2, 2), (1, 2)]])
    with pytest.raises(ValueError, match="holes"):
        to_polygon(ring)


def test_empty_polygon_is_rejected():
    with pytest.raises(ValueError, match="empty"):
        to_polygon(Polygon())


def test_to_dict_is_json_ready():
    d = describe(HEXAGON).to_dict()
    assert set(d) == {"area", "perimeter", "nVertices", "inradius", "incenter", "circumradius",
                      "length", "width", "elongation", "axisAngle", "convexity"}
    assert json.loads(json.dumps(d)) == d
