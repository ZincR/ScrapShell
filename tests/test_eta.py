"""Tests for sizing/eta.py: per-face eta files for PH-CPF runs."""
import numpy as np
import pytest

from sizing.eta import (ETA_MAX, eta_for_area, face_adjacency, limit_aspect, log_jump, main,
                        smooth_eta, smooth_log, write_eta)
from sizing.footprint import DEFAULT_ETA, footprint, load_surface, panel_count


@pytest.fixture(scope="module")
def metro():
    return load_surface("metro")


@pytest.mark.parametrize("max_aspect", [1.5, 2.0, 3.0])
def test_limit_aspect(metro, max_aspect):
    eta = limit_aspect(metro, 4 * DEFAULT_ETA, max_aspect)
    free = footprint(metro, eta=4 * DEFAULT_ETA)
    fp = footprint(metro, eta=eta)
    assert eta.shape == (metro.n_faces,) and np.all(eta <= 4 * DEFAULT_ETA)
    assert fp.aspect.max() <= max_aspect * (1 + 1e-9)
    over = free.aspect > max_aspect
    assert over.any() and np.allclose(fp.aspect[over], max_aspect)     # lowered just enough
    assert np.array_equal(eta[~over], np.full((~over).sum(), 4 * DEFAULT_ETA))


def test_limit_aspect_per_face_and_bad_input(metro):
    base = np.full(metro.n_faces, DEFAULT_ETA)
    assert np.array_equal(limit_aspect(metro, base, 100.0), base)
    with pytest.raises(ValueError, match="greater than 1"):
        limit_aspect(metro, DEFAULT_ETA, 1.0)
    with pytest.raises(ValueError, match="one per face"):
        limit_aspect(metro, np.full(3, DEFAULT_ETA), 2.0)


def test_write_eta_and_command_line(metro, tmp_path, capsys):
    path = write_eta(tmp_path / "sub" / "eta.txt", np.array([0.002, 0.0031415926535]))
    assert np.allclose(np.loadtxt(path), [0.002, 0.0031415926535], rtol=1e-9)
    out = tmp_path / "metro.txt"
    assert main(["metro", "--eta", "0.008", "--max-aspect", "2", "--out", str(out)]) == 0
    eta = np.loadtxt(out)
    assert np.allclose(eta, limit_aspect(metro, 0.008, 2.0), rtol=1e-9)
    assert "below 0.008 on" in capsys.readouterr().out


def test_face_adjacency(metro):
    adj = face_adjacency(metro.faces)
    assert adj.shape == (2087, 2) and np.all(adj[:, 0] != adj[:, 1])
    shared = [len(set(metro.faces[f]) & set(metro.faces[g])) for f, g in adj]
    assert set(shared) == {2}
    assert np.array_equal(face_adjacency(np.array([[0, 1, 2], [2, 1, 3], [3, 4, 5]])), [[0, 1]])


def test_log_jump_and_smooth_log():
    adj = np.array([[0, 1], [1, 2]])
    assert np.allclose(log_jump([1.0, np.e, np.e], adj), [1.0, 0.0])
    v = np.array([1.0, 8.0, 1.0])
    assert np.allclose(smooth_log(v, adj, 0), v)
    once = smooth_log(v, adj, 1, lam=1.0)                       # each face takes its neighbours' mean
    assert np.allclose(np.log(once), [np.log(8), 0.0, np.log(8)])
    many = smooth_log(v, adj, 500)
    assert np.ptp(np.log(many)) < 1e-6                        # a connected graph smooths to a constant
    with pytest.raises(ValueError):
        smooth_log(v, adj, -1)


def test_eta_for_area_inverts_the_footprint(metro):
    rng = np.random.default_rng(0)
    eta = DEFAULT_ETA * np.exp(rng.uniform(np.log(0.5), np.log(8), metro.n_faces))
    target = footprint(metro, eta=eta).panel_area()
    assert np.allclose(eta_for_area(metro, target), eta, rtol=1e-6)
    # a target the curvature cannot reach stops at eta_max, and the panel is as large as it gets
    big = eta_for_area(metro, 1e6)
    assert np.allclose(big, ETA_MAX)
    with pytest.raises(ValueError):
        eta_for_area(metro, -1.0)


def test_eta_for_area_with_max_aspect(metro):
    target = footprint(metro, eta=4 * DEFAULT_ETA).panel_area()
    eta = eta_for_area(metro, target, max_aspect=2.0)
    assert np.allclose(eta, limit_aspect(metro, 4 * DEFAULT_ETA, 2.0), rtol=1e-6)
    assert footprint(metro, eta=eta).aspect.max() <= 2 * (1 + 1e-9)


def test_smooth_eta(metro):
    adj = face_adjacency(metro.faces)
    assert np.allclose(smooth_eta(metro, 4 * DEFAULT_ETA, steps=0), 4 * DEFAULT_ETA, rtol=1e-6)
    rough = footprint(metro, eta=4 * DEFAULT_ETA)
    fp = footprint(metro, eta=smooth_eta(metro, 4 * DEFAULT_ETA, steps=400))
    size = lambda f: np.sqrt(f.a1 * f.a2)
    assert np.median(log_jump(size(fp), adj)) < 0.3 * np.median(log_jump(size(rough), adj))
    assert abs(panel_count(metro, fp) / panel_count(metro, rough) - 1) < 0.1   # about as many panels
    capped = footprint(metro, eta=smooth_eta(metro, 4 * DEFAULT_ETA, steps=400, max_aspect=2.0))
    assert capped.aspect.max() <= 2 * (1 + 1e-9)


def test_command_line_smooth(metro, tmp_path, capsys):
    out = tmp_path / "smooth.txt"
    assert main(["metro", "--eta", "0.008", "--smooth", "50", "--out", str(out)]) == 0
    assert np.allclose(np.loadtxt(out), smooth_eta(metro, 0.008, 50), rtol=1e-9)
    assert "log-size jump" in capsys.readouterr().out
