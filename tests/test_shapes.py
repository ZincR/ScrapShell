"""Tests for stock/shapes.py: every family yields usable scraps in the canonical frame."""
import math
from collections import Counter

import numpy as np
import pytest
from shapely.geometry import Point, box

from stock.describe import INRADIUS_REL_TOL, describe
from stock.shapes import (FAMILIES, Limits, blob, check, generate_shape, interior_angles,
                          sample_shapes)

N_PER_FAMILY = 1000


@pytest.fixture(scope="module", params=sorted(FAMILIES))
def batch(request):
    """N_PER_FAMILY scraps of one family, their descriptors, and the draw outcomes."""
    family = request.param
    rng = np.random.default_rng(sorted(FAMILIES).index(family))
    stats = Counter()
    shapes = [generate_shape(rng, family, stats=stats) for _ in range(N_PER_FAMILY)]
    return family, shapes, [describe(p) for p in shapes], stats


# --- every family -----------------------------------------------------------------------------

def test_all_scraps_pass_the_checks(batch):
    family, shapes, _, _ = batch
    reasons = Counter(check(p) for p in shapes)
    assert reasons == {None: N_PER_FAMILY}, f"{family}: {reasons}"


def test_canonical_frame(batch):
    """Unit area, an inscribed circle centred at the origin, long axis along x."""
    family, shapes, descs, _ = batch
    for p, d in zip(shapes, descs):
        assert d.area == pytest.approx(1.0, rel=1e-9)
        assert p.contains(Point(0, 0))
        # the origin is a centre of a largest inscribed circle (two estimates, each within tolerance)
        assert p.exterior.distance(Point(0, 0)) == pytest.approx(d.inradius, abs=2.5 * INRADIUS_REL_TOL)
        # the axis-aligned bounding box is a minimum-area enclosing rectangle, wider than tall
        # (this holds even when several rectangles tie for the minimum, as for triangles)
        x0, y0, x1, y1 = p.bounds
        assert (x1 - x0) * (y1 - y0) == pytest.approx(d.length * d.width, rel=1e-9), family
        assert x1 - x0 >= y1 - y0 - 1e-9, family


def test_interior_angles_sum(batch):
    _, shapes, _, _ = batch
    for p in shapes:
        pts = np.asarray(p.exterior.coords)[:-1]
        assert interior_angles(pts).sum() == pytest.approx((len(pts) - 2) * 180.0)


def test_acceptance_rate(batch):
    """Most draws pass the checks, so generating a stock stays cheap."""
    family, _, _, stats = batch
    assert stats["ok"] == N_PER_FAMILY
    assert stats["ok"] / sum(stats.values()) > 0.5, f"{family}: {dict(stats)}"


# --- family character -------------------------------------------------------------------------

def test_family_character(batch):
    family, shapes, descs, _ = batch
    convexity = np.array([d.convexity for d in descs])
    if family == "board":
        for p, d in zip(shapes, descs):
            assert d.n_vertices == 4
            assert interior_angles(np.asarray(p.exterior.coords)[:-1]) == pytest.approx([90] * 4)
            assert 1.0 <= d.elongation <= 5.0 + 1e-9
    elif family == "shard":
        assert (convexity > 1 - 1e-9).all()
    elif family == "remnant":
        assert np.mean(convexity < 0.99) > 0.95
    elif family == "blob":
        assert np.mean(convexity < 0.99) > 0.3


# --- check() ----------------------------------------------------------------------------------

L_SHAPE = [(0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2)]


@pytest.mark.parametrize("outline, kind", [
    (L_SHAPE, None),
    (box(0, 0, 10, 1), "fatness"),                                            # a 10:1 strip
    (box(0, 0, 1, 1).union(box(1, 0.48, 2, 0.52)).union(box(2, 0, 3, 1)), "feature"),   # dumbbell
    ([(0, 0), (0.45, 0), (0.45, 0.6), (0.55, 0.6), (0.55, 0), (1, 0), (1, 1), (0, 1)], "notch"),  # slot
    ([(0, 0), (1, 0), (0.5, 0.5 / math.tan(math.radians(5)))], "angle"),      # 10 degree spike
    ([(0, 0), (1, 0), (1, 1), (0.01, 1), (0, 0.99)], "edge"),                 # chipped corner
    ([(math.cos(a), math.sin(a)) for a in np.linspace(0, 2 * math.pi, 20, endpoint=False)], "vertices"),
    ([(0, 0), (1, 1), (1, 0), (0, 1)], "outline"),                            # bow-tie
])
def test_check_reasons(outline, kind):
    reason = check(outline)
    assert (reason if reason is None else reason.split(":")[0]) == kind


def test_limits_are_adjustable():
    twenty_gon = [(math.cos(a), math.sin(a)) for a in np.linspace(0, 2 * math.pi, 20, endpoint=False)]
    assert check(twenty_gon, Limits(max_vertices=24)) is None
    assert check(box(0, 0, 10, 1), Limits(min_fatness=0.1)) is None


# --- sampling ---------------------------------------------------------------------------------

def test_reproducible():
    a = sample_shapes(np.random.default_rng(7), 30)
    b = sample_shapes(np.random.default_rng(7), 30)
    c = sample_shapes(np.random.default_rng(8), 30)
    assert [f for f, _ in a] == [f for f, _ in b]
    assert all(np.array_equal(p.exterior.coords, q.exterior.coords) for (_, p), (_, q) in zip(a, b))
    assert not all(p.equals(q) for (_, p), (_, q) in zip(a, c))


def test_mix_frequencies():
    picks = Counter(f for f, _ in sample_shapes(np.random.default_rng(0), 2000, {"board": 3, "blob": 1}))
    sigma = math.sqrt(2000 * 0.75 * 0.25)
    assert abs(picks["board"] - 1500) < 4 * sigma
    assert picks["board"] + picks["blob"] == 2000


@pytest.mark.parametrize("mix", [{"board": -1, "blob": 1}, {"board": 0}, {"plank": 1}])
def test_bad_mix(mix):
    with pytest.raises(ValueError):
        sample_shapes(np.random.default_rng(0), 5, mix)


def test_unknown_family():
    with pytest.raises(ValueError, match="unknown family"):
        generate_shape(np.random.default_rng(0), "plank")


def test_impossible_limits_give_up():
    with pytest.raises(RuntimeError, match="no usable"):
        generate_shape(np.random.default_rng(0), "blob", limits=Limits(min_fatness=0.9), max_tries=20)


def test_blob_irregularity_bound():
    with pytest.raises(ValueError, match="irregularity"):
        blob(np.random.default_rng(0), irregularity=1.0)


def test_family_params_pass_through():
    rng = np.random.default_rng(0)
    for _ in range(50):
        d = describe(generate_shape(rng, "board", {"aspect_median": 1.2, "aspect_max": 1.5}))
        assert d.elongation <= 1.5 + 1e-9
