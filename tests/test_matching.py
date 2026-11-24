"""Tests for matching/: the fit test and the assignment."""
import numpy as np
import pytest
from shapely.geometry import Polygon

from matching.assign import assign, max_matching
from matching.fit import (centered, describe_panel, fit, fit_graph, load_deposit, min_width,
                          model_panels, rotate)


def rect(w, h):
    return np.array([[0, 0], [w, 0], [w, h], [0, h]], float)


def test_min_width():
    assert min_width(rect(2, 0.5)) == pytest.approx(0.5)
    tri = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]])
    assert min_width(tri) == pytest.approx(np.sqrt(3) / 2)


def test_fit_translation_and_rotation():
    assert fit(Polygon(rect(1.05, 1.05)), centered(rect(1.0, 0.2))) is not None
    assert fit(Polygon(rect(0.95, 0.95)), centered(rect(1.0, 0.2)), rotations=1) is None
    r = fit(Polygon(rect(1, 1)), centered(rect(1.2, 0.1)))      # only along a diagonal
    assert r is not None and 0 < r[0] % (np.pi / 2) < np.pi / 2


def test_fit_non_convex_scrap_and_panel():
    L = Polygon([(0, 0), (2, 0), (2, 0.6), (0.6, 0.6), (0.6, 2), (0, 2)])
    assert fit(L, centered(rect(0.5, 0.5))) is not None
    assert fit(L, centered(rect(0.7, 0.7))) is None
    notched = Polygon([(0, 0), (1, 0), (1, 1), (0.55, 1), (0.5, 0.2), (0.45, 1), (0, 1)])
    assert fit(notched, centered(rect(0.8, 0.8))) is None      # corners inside, notch through it
    bowtie = np.array([[0, 0], [0.6, 0.15], [1.2, 0], [1.2, 0.6], [0.6, 0.45], [0, 0.6]])
    S = Polygon(rect(1.25, 0.65))
    angle, t = fit(S, centered(bowtie))
    assert S.buffer(1e-9).contains(Polygon(rotate(centered(bowtie), angle) + t))
    assert fit(Polygon(rect(1.25, 0.5)), centered(bowtie)) is None


def test_assign_most_panels_then_least_waste():
    fits = np.array([[1, 1, 0],
                     [1, 0, 0],
                     [0, 0, 0]], bool)
    a = assign(fits, np.array([1.0, 1.0, 1.0]), np.array([2.0, 5.0, 1.5]))
    assert list(a.scrap) == [1, 0, -1]                           # panel 0 yields scrap 0 to panel 1
    assert a.fit_rate == pytest.approx(2 / 3) and a.utilization == pytest.approx(2 / 7)
    fits = np.array([[1, 1], [1, 1]], bool)
    a = assign(fits, np.array([1.0, 2.0]), np.array([3.0, 2.5]))
    assert set(a.scrap) == {0, 1} and a.assigned.sum() == max_matching(fits) == 2


def test_metro_panels_and_a_small_graph():
    panels = model_panels("metro")
    assert len(panels) == 436
    assert np.median([Polygon(P).area for P in panels]) == pytest.approx(0.258, abs=0.005)
    deposit = load_deposit()
    small = {"scraps": deposit["scraps"][:300]}
    g = fit_graph(panels[:5], small, rotations=12, workers=1)
    for i, pl in enumerate(g.placements):
        P = centered(panels[i])
        for j, (angle, t) in pl.items():                      # every placement is really inside
            S = Polygon(small["scraps"][j]["polygon"])
            assert S.buffer(1e-6).contains(Polygon(rotate(P, angle) + t))
        d = describe_panel(panels[i])
        for j in range(300):                                  # sure fits are always found
            if small["scraps"][j]["inradius"] >= d.circumradius:
                assert j in pl
