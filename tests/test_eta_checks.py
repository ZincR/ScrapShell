"""The per-face eta hook in PH-CPF, checked on the committed Metro runs of matlab/run_eta_checks.m."""
import json

import numpy as np
import pytest

from sizing.footprint import DEFAULT_ETA, footprint, load_surface, run_eta
from stock.model_ref import WEB_DATA, decode_array
from stock.place import model_ids

CHECKS = ["metro_eta_1x", "metro_eta_4x", "metro_eta_ramp", "metro_eta_4x_a2", "metro_eta_4x_delta",
          "metro_eta_4x_smooth", "metro_eta_4x_smooth_a2"]
PER_FACE = ["metro_eta_1x", "metro_eta_4x", "metro_eta_ramp", "metro_eta_4x_a2", "metro_eta_4x_smooth",
            "metro_eta_4x_smooth_a2"]


def _data(run):
    return json.loads((WEB_DATA / f"{run}.json").read_text())


def _meta(run):
    return json.loads((WEB_DATA / f"{run}.meta.json").read_text())


def test_variants_stay_out_of_the_model_list():
    assert not set(CHECKS) & set(model_ids())
    assert set(CHECKS) <= set(model_ids(variants=True))
    assert len(model_ids()) == 8


@pytest.mark.parametrize("run", CHECKS)
def test_runs_are_metro_variants(run):
    d, m = _data(run), _meta(run)
    assert d["id"] == m["id"] == run and d["model"] == m["model"] == "metro"
    assert m["group"] == "Sizing checks"
    eta = run_eta(d)
    if run in PER_FACE:
        assert isinstance(eta, np.ndarray) and eta.shape == (m["nf"],)
    else:
        assert eta == pytest.approx(4 * DEFAULT_ETA) and d["params"]["delta"] == pytest.approx(0.02)


def test_outcomes():
    """What docs/plan.md and docs/progress.md report: the stage each run reached, and the
    continuation level (iter_opt) of its last attempt. Rounding fails at every level for 4x eta,
    with or without the aspect cap; it succeeds at once when delta grows with eta, and after some
    doublings when the 4x size field is smoothed."""
    expected = {"metro_eta_1x": (8, 1), "metro_eta_4x": (4, 10), "metro_eta_ramp": (8, 5),
                "metro_eta_4x_a2": (4, 10), "metro_eta_4x_delta": (8, 0),
                "metro_eta_4x_smooth": (8, 4), "metro_eta_4x_smooth_a2": (8, 8)}
    for run, (stage, level) in expected.items():
        assert (max(_meta(run)["stages"]), _data(run)["params"]["iter_opt"]) == (stage, level), run


@pytest.mark.parametrize("run", CHECKS)
def test_hook_matches_the_footprint_model(run):
    """MATLAB's per-face sizing ellipses equal the Python model's, from the baseline's curvatures."""
    surf = load_surface("metro")
    d = _data(run)
    fp = footprint(surf, d["params"]["delta"], run_eta(d))
    st = d["stages"]["guide"]
    for a, key in ((fp.a1, "axis1"), (fp.a2, "axis2")):
        stored = np.linalg.norm(decode_array(st[key]), axis=1) * surf.meters_per_unit
        assert np.allclose(a, stored, rtol=2e-5)


def test_the_ramp_spans_1x_to_16x_along_the_longest_axis():
    surf = load_surface("metro")
    eta = run_eta(_data("metro_eta_ramp"))
    # set from vertex extremes, sampled at face centroids: just inside [1x, 16x]
    assert DEFAULT_ETA * (1 - 1e-6) <= eta.min() < 1.1 * DEFAULT_ETA
    assert 14 * DEFAULT_ETA < eta.max() <= 16 * DEFAULT_ETA * (1 + 1e-6)
    ax = np.argmax(np.ptp(surf.vertices, axis=0))
    order = np.argsort(surf.centroids[:, ax])
    assert np.all(np.diff(np.log(eta[order])) >= -1e-6)


def test_eta_per_face_equal_to_the_baseline_reproduces_it():
    base, same = _data("metro"), _data("metro_eta_1x")
    assert np.allclose(run_eta(same), DEFAULT_ETA)
    for stage, key in (("guide", "axis1"), ("dual", "idx"), ("dual", "D")):
        assert np.array_equal(decode_array(same["stages"][stage][key]),
                              decode_array(base["stages"][stage][key]))
    Vb = decode_array(base["stages"]["planar"]["frames"][-1]["V"])
    Vs = decode_array(same["stages"]["planar"]["frames"][-1]["V"])
    assert np.array_equal(Vb, Vs)
    mb, ms = _meta("metro")["metrics"], _meta("metro_eta_1x")["metrics"]
    for k in ("nSing", "nPanels", "planMax", "distMax"):
        assert ms[k] == mb[k]


@pytest.mark.parametrize("run,args", [("metro_eta_4x_smooth", {}),
                                      ("metro_eta_4x_smooth_a2", {"max_aspect": 2.0})])
def test_smoothed_runs_use_the_smoothed_size_field(run, args):
    """The eta of each smoothed run is what sizing.eta.smooth_eta gives (400 steps from 4x), and its
    footprints are smoother than those of plain 4x eta."""
    from sizing.eta import face_adjacency, log_jump, smooth_eta
    surf = load_surface("metro")
    eta = run_eta(_data(run))
    assert np.allclose(eta, smooth_eta(surf, 4 * DEFAULT_ETA, 400, **args), rtol=1e-6)
    adj = face_adjacency(surf.faces)
    jump = lambda e: np.median(log_jump(np.sqrt(footprint(surf, eta=e).a1 * footprint(surf, eta=e).a2), adj))
    assert jump(eta) < 0.5 * jump(4 * DEFAULT_ETA)
