"""Tests for stock/model_ref.py, against the committed viewer data in web/data/."""
import base64

import numpy as np
import pytest
import shapely

from stock.model_ref import decode_array, flatten_panel, load_model_ref


@pytest.fixture(scope="module")
def metro():
    return load_model_ref("metro")


def test_metro_reference(metro):
    """Matches the numbers measured from the baseline run (docs/progress.md)."""
    assert metro.surface_area == pytest.approx(178.276, rel=1e-4)
    assert metro.bbd == pytest.approx(17.308, rel=1e-4)
    assert len(metro.panels) == 436
    assert metro.a_ref == pytest.approx(0.34286, rel=1e-4)
    assert np.median(metro.panel_circumradii) == pytest.approx(0.4446, rel=1e-3)
    # an open shell: the panels cover the surface minus the stripped boundary layer
    assert 0.8 < metro.panel_areas.sum() / metro.surface_area < 0.9


def test_to_dict(metro):
    d = metro.to_dict()
    assert d["source"] == "web/data/metro.json"
    assert d["nPanels"] == 436
    assert d["panelAreaMedian"] == metro.a_ref


def test_panels_are_flat_and_counter_clockwise(metro):
    for p in metro.panels:
        assert p.shape[1] == 2 and len(p) >= 3
        assert shapely.is_ccw(shapely.LinearRing(p))
    assert np.mean([len(p) == 6 for p in metro.panels]) > 0.9      # mostly hexagons (422 of 436)


def test_flatten_panel_keeps_area_and_orients():
    rng = np.random.default_rng(0)
    hexagon = np.c_[np.cos(np.arange(6) * np.pi / 3), np.sin(np.arange(6) * np.pi / 3), np.zeros(6)]
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))          # random rotation (or reflection)
    for pts in (hexagon @ Q.T + [1, 2, 3], (hexagon @ Q.T)[::-1]):
        q = flatten_panel(pts)
        assert shapely.Polygon(q).area == pytest.approx(3 * np.sqrt(3) / 2)
        assert shapely.is_ccw(shapely.LinearRing(q))


def test_decode_array_round_trip():
    x = np.arange(12, dtype=np.float32).reshape(4, 3) / 7
    packed = {"dtype": "f32", "shape": [4, 3], "b64": base64.b64encode(x.tobytes()).decode()}
    assert np.array_equal(decode_array(packed), x.astype(np.float64))
    i = np.array([[0, 5], [7, 2]], dtype=np.int32)
    packed = {"dtype": "i32", "shape": [2, 2], "b64": base64.b64encode(i.tobytes()).decode()}
    assert np.array_equal(decode_array(packed), i)


def test_model_without_panels():
    """car_r stopped before planarization: its surface loads, its panels are empty."""
    car = load_model_ref("car_r")
    assert car.panels == () and car.surface_area > 0 and car.bbd > 0
    assert "panelAreaMedian" not in car.to_dict()
    with pytest.raises(ValueError, match="no planarized panels"):
        car.a_ref


def test_unknown_model():
    with pytest.raises(ValueError, match="no viewer data"):
        load_model_ref("teapot")
