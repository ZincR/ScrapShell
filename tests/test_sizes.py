"""Tests for stock/sizes.py: samples follow the requested (truncated) distributions."""
import math

import numpy as np
import pytest
from scipy import stats

from stock import sizes

UNIFORM = {"kind": "uniform", "min": 1.0, "max": 6.0}
LOGNORMAL = {"kind": "lognormal", "median": 2.5, "sigma": 0.5, "min": 0.5, "max": 12.0}
LOGNORMAL_TIGHT = {"kind": "lognormal", "median": 2.5, "sigma": 0.8, "min": 2.0, "max": 3.0}
PARETO = {"kind": "pareto", "alpha": 1.5, "min": 1.0, "max": 20.0}
MIXTURE = {"kind": "mixture", "components": [
    {"weight": 0.3, "dist": {"kind": "uniform", "min": 1.0, "max": 2.0}},
    {"weight": 0.7, "dist": {"kind": "uniform", "min": 5.0, "max": 6.0}},
]}


def truncated_cdf(dist):
    lo, hi = dist["min"], dist["max"]
    if dist["kind"] == "uniform":
        return stats.uniform(lo, hi - lo).cdf
    if dist["kind"] == "pareto":
        a = dist["alpha"]
        return lambda x: (1 - (lo / np.asarray(x)) ** a) / (1 - (lo / hi) ** a)
    mu, s = math.log(dist["median"]), dist["sigma"]
    Fa, Fb = stats.norm.cdf((math.log(lo) - mu) / s), stats.norm.cdf((math.log(hi) - mu) / s)
    return lambda x: (stats.norm.cdf((np.log(x) - mu) / s) - Fa) / (Fb - Fa)


@pytest.mark.parametrize("dist", [UNIFORM, LOGNORMAL, LOGNORMAL_TIGHT, PARETO],
                         ids=["uniform", "lognormal", "lognormal-tight", "pareto"])
def test_samples_follow_the_distribution(dist):
    x = sizes.sample(np.random.default_rng(0), dist, 20000)
    assert len(x) == 20000
    assert x.min() >= dist["min"] and x.max() <= dist["max"]
    assert stats.kstest(x, truncated_cdf(dist)).pvalue > 1e-3


def test_mixture_weights():
    x = sizes.sample(np.random.default_rng(0), MIXTURE, 20000)
    small = np.mean(x < 3)
    assert abs(small - 0.3) < 4 * math.sqrt(0.3 * 0.7 / 20000)
    assert ((1 <= x) & (x <= 2) | (5 <= x) & (x <= 6)).all()


def test_reproducible():
    a = sizes.sample(np.random.default_rng(5), MIXTURE, 100)
    b = sizes.sample(np.random.default_rng(5), MIXTURE, 100)
    assert np.array_equal(a, b)


def test_zero_samples_consume_no_randomness():
    rng = np.random.default_rng(1)
    assert len(sizes.sample(rng, LOGNORMAL, 0)) == 0
    assert rng.random() == np.random.default_rng(1).random()


@pytest.mark.parametrize("dist, message", [
    ({"kind": "gamma"}, "unknown"),
    ({"kind": "uniform", "min": 3, "max": 2}, "min < max"),
    ({"kind": "uniform", "min": 0, "max": 2}, "min < max"),
    ({"kind": "uniform", "min": 1}, "min < max"),
    ({"kind": "lognormal", "median": 2, "sigma": 0, "min": 1, "max": 3}, "sigma"),
    ({"kind": "lognormal", "sigma": 0.5, "min": 1, "max": 3}, "median"),
    ({"kind": "pareto", "alpha": 0, "min": 1, "max": 3}, "alpha"),
    ({"kind": "mixture", "components": []}, "at least one"),
    ({"kind": "mixture", "components": [{"weight": -1, "dist": UNIFORM}]}, "weights"),
    ({"kind": "mixture", "components": [{"weight": 0, "dist": UNIFORM}]}, "weights"),
    ({"kind": "mixture", "components": [{"weight": 1, "dist": {"kind": "gamma"}}]}, "unknown"),
])
def test_malformed_distributions(dist, message):
    with pytest.raises(ValueError, match=message):
        sizes.sample(np.random.default_rng(0), dist, 10)


def test_lognormal_outside_its_bounds():
    dist = {"kind": "lognormal", "median": 1.0, "sigma": 0.1, "min": 100.0, "max": 101.0}
    with pytest.raises(ValueError, match="widen the bounds"):
        sizes.sample(np.random.default_rng(0), dist, 10)
