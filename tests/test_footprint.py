"""Tests for sizing/footprint.py, against the committed viewer data in web/data/."""
import json
import math

import numpy as np
import pytest

from sizing.footprint import (DEFAULT_DELTA, DEFAULT_ETA, HEX_AREA, PANEL_K, Footprint, calibrate,
                              capped_fraction, footprint, load_surface, main, panel_count)
from stock.model_ref import WEB_DATA, decode_array
from stock.place import model_ids, placement

COMPLETED = ["metro", "intersection", "six_s0", "vase", "sphere2con_s0", "c_s0"]


@pytest.fixture(scope="module")
def metro():
    return load_surface("metro")


def test_surface_in_meters(metro):
    pl = placement("metro")
    assert metro.span == 15.0 and metro.meters_per_unit == pytest.approx(pl.meters_per_unit)
    assert metro.area == pytest.approx(pl.surface_area)
    assert np.linalg.norm(metro.vertices.max(0) - metro.vertices.min(0)) == pytest.approx(15.0)
    assert metro.n_faces == len(metro.face_area) == len(metro.kmin) == 1418


def test_principal_directions(metro):
    normals = np.cross(*(metro.vertices[metro.faces[:, i]] - metro.vertices[metro.faces[:, 0]]
                         for i in (1, 2)))
    normals /= np.linalg.norm(normals, axis=1, keepdims=True)
    for d in (metro.dmin, metro.dmax):
        assert np.allclose(np.linalg.norm(d, axis=1), 1, atol=1e-6)
        assert np.abs(np.sum(d * normals, axis=1)).max() < 1e-3
    assert np.abs(np.sum(metro.dmin * metro.dmax, axis=1)).max() < 1e-3
    assert (metro.kmin >= 0).all() and (metro.kmax >= 0).all()


@pytest.mark.parametrize("model", model_ids())
def test_reproduces_the_baseline_ellipses(model):
    """At the baseline's delta and eta, the recovered curvatures give back the stored semi-axes."""
    surf = load_surface(model)
    fp = footprint(surf)
    st = json.loads((WEB_DATA / f"{model}.json").read_text())["stages"]["guide"]
    for a, key in ((fp.a1, "axis1"), (fp.a2, "axis2")):
        stored = np.linalg.norm(decode_array(st[key]), axis=1) * surf.meters_per_unit
        assert np.allclose(a, stored, rtol=1e-5)


def test_scale_invariance():
    small, big = load_surface("vase", span=15.0), load_surface("vase", span=30.0)
    fs, fb = footprint(small), footprint(big)
    assert np.allclose(fb.a1, 2 * fs.a1) and np.allclose(fb.a2, 2 * fs.a2)
    assert np.allclose(big.kmin, small.kmin / 2)
    assert panel_count(big, fb) == pytest.approx(panel_count(small, fs))
    assert np.array_equal(fb.capped, fs.capped)


def test_eta_caps_the_footprint(metro):
    fp = footprint(metro)
    assert fp.delta == pytest.approx(DEFAULT_DELTA * 15.0)
    assert fp.eta == pytest.approx(DEFAULT_ETA * metro.area)
    assert fp.major.max() <= math.sqrt(fp.eta) * (1 + 1e-12)
    # without the cap, the semi-axes are sqrt(delta / |k|)
    free = footprint(metro, eta=1e6)
    curved = metro.kmin > 1e-3
    assert np.allclose(free.a1[curved], np.sqrt(free.delta / metro.kmin[curved]), rtol=1e-6)


def test_larger_eta_means_fewer_larger_panels(metro):
    fps = [footprint(metro, eta=e) for e in (0.002, 0.008, 0.032)]
    counts = [panel_count(metro, fp) for fp in fps]
    capped = [capped_fraction(metro, fp) for fp in fps]
    assert counts[0] > counts[1] > counts[2]
    assert capped[0] > capped[1] > capped[2]
    assert all(np.all(b.major >= a.major) for a, b in zip(fps, fps[1:]))
    # ... and more stretched: the aspect rises toward sqrt(k_max / k_min)
    assert all(np.all(b.aspect >= a.aspect * (1 - 1e-12)) for a, b in zip(fps, fps[1:]))
    ks, kb = np.minimum(metro.kmin, metro.kmax), np.maximum(metro.kmin, metro.kmax)
    curved = ks > 1e-3
    assert np.all(fps[-1].aspect[curved] <= np.sqrt(kb[curved] / ks[curved]) * (1 + 1e-9))


def test_panel_area():
    fp = Footprint(a1=np.array([1.0, 2.0]), a2=np.array([1.0, 0.5]), delta=0.1, eta=4.0,
                   capped=np.array([True, False]))
    assert np.allclose(fp.panel_area(), HEX_AREA * PANEL_K**2)
    assert np.allclose(fp.panel_radius(), PANEL_K * np.array([1.0, 2.0]))
    assert np.allclose(fp.aspect, [1.0, 4.0])


def test_eta_per_face(metro):
    """One eta per face gives each face the footprint it would get from that eta alone."""
    base = footprint(metro)
    same = footprint(metro, eta=np.full(metro.n_faces, DEFAULT_ETA))
    assert np.array_equal(same.a1, base.a1) and np.array_equal(same.a2, base.a2)
    assert np.array_equal(same.capped, base.capped)

    eta = np.where(metro.region == 0, 0.032, 0.002)          # larger panels where the surface is flat
    fp = footprint(metro, eta=eta)
    big = footprint(metro, eta=0.032)
    flat = metro.region == 0
    assert flat.any() and (~flat).any()
    assert np.allclose(fp.a1[flat], big.a1[flat]) and np.allclose(fp.a2[~flat], base.a2[~flat])
    assert fp.eta.shape == (metro.n_faces,)
    assert np.allclose(fp.eta, eta * metro.area)
    assert panel_count(metro, big) < panel_count(metro, fp) < panel_count(metro, base)
    # never larger than the curvature allows
    curved = np.minimum(metro.kmin, metro.kmax) > 1e-3
    assert np.all(fp.a1[curved] <= np.sqrt(fp.delta / metro.kmin[curved]) * (1 + 1e-9))


@pytest.mark.parametrize("model", COMPLETED)
def test_panel_k_matches_the_baseline(model):
    """One constant predicts every model's panels: their size, and how many there are."""
    c = calibrate(model)
    assert c.k == pytest.approx(PANEL_K, rel=0.05)
    assert c.k_p10 < c.k < c.k_p90
    assert c.predicted == pytest.approx(c.n_panels, rel=0.08)
    assert 0.7 < c.covered <= 1.0


@pytest.mark.parametrize("model", ["car_r", "3holes_r"])
def test_models_without_panels(model):
    surf = load_surface(model)
    assert panel_count(surf, footprint(surf)) > 0
    assert calibrate(model) is None


def test_bad_input():
    with pytest.raises(ValueError, match="positive"):
        load_surface("metro", span=0)
    vase = load_surface("vase")
    with pytest.raises(ValueError, match="positive"):
        footprint(vase, eta=0)
    with pytest.raises(ValueError, match="positive"):
        footprint(vase, eta=np.r_[0.002 * np.ones(vase.n_faces - 1), -1.0])
    with pytest.raises(ValueError, match="one per face"):
        footprint(vase, eta=np.full(vase.n_faces + 1, 0.002))
    with pytest.raises(ValueError, match="no viewer data"):
        load_surface("teapot")


def test_command_line(capsys):
    assert main(["--models", "metro", "car_r"]) == 0
    out = capsys.readouterr().out
    assert "metro" in out and "car_r" in out and "no baseline panels" in out
