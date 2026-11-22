"""Distributions of a positive scrap measure. The stock draws scrap size (sqrt(area), in meters).

A distribution is a plain dict, so it can be stored in JSON with the stock:

    {"kind": "uniform",   "min": 0.2, "max": 1.0}
    {"kind": "lognormal", "median": 0.5, "sigma": 0.3, "min": 0.15, "max": 1.5}
    {"kind": "pareto",    "alpha": 1.5, "min": 0.2, "max": 1.5}       # heavy tail, truncated
    {"kind": "mixture",   "components": [{"weight": 0.5, "dist": {...}}, ...]}

min and max bound every kind. A lognormal is conditioned on [min, max] by redrawing, so there is
no pile-up at the bounds. The truncated Pareto is sampled exactly by inverting its CDF.
"""
from __future__ import annotations

import math
from collections.abc import Mapping

import numpy as np

KINDS = ("uniform", "lognormal", "pareto", "mixture")


def validate(dist: Mapping) -> None:
    """Raise ValueError if dist is not a well-formed size distribution."""
    kind = dist.get("kind")
    if kind not in KINDS:
        raise ValueError(f"unknown size distribution {kind!r}; choose from {list(KINDS)}")
    if kind == "mixture":
        comps = dist.get("components") or []
        if not comps:
            raise ValueError("a mixture needs at least one component")
        weights = [c.get("weight", -1) for c in comps]
        if min(weights) < 0 or sum(weights) <= 0:
            raise ValueError("mixture weights must be non-negative and not all zero")
        for c in comps:
            validate(c["dist"])
        return
    lo, hi = dist.get("min"), dist.get("max")
    if lo is None or hi is None or not 0 < lo < hi:
        raise ValueError(f"{kind}: need 0 < min < max, got min={lo}, max={hi}")
    if kind == "lognormal" and not (dist.get("median", 0) > 0 and dist.get("sigma", 0) > 0):
        raise ValueError("lognormal: need median > 0 and sigma > 0")
    if kind == "pareto" and not dist.get("alpha", 0) > 0:
        raise ValueError("pareto: need alpha > 0")


def sample(rng: np.random.Generator, dist: Mapping, n: int) -> np.ndarray:
    """n values drawn from dist."""
    validate(dist)
    if n == 0:
        return np.empty(0)
    kind = dist["kind"]
    if kind == "mixture":
        comps = dist["components"]
        w = np.array([c["weight"] for c in comps], dtype=float)
        which = rng.choice(len(comps), size=n, p=w / w.sum())
        out = np.empty(n)
        for k, c in enumerate(comps):
            mask = which == k
            out[mask] = sample(rng, c["dist"], int(mask.sum()))
        return out
    lo, hi = float(dist["min"]), float(dist["max"])
    if kind == "uniform":
        return rng.uniform(lo, hi, n)
    if kind == "pareto":
        a = float(dist["alpha"])
        u = rng.uniform(0.0, 1.0, n)
        return lo * (1.0 - u * (1.0 - (lo / hi) ** a)) ** (-1.0 / a)
    # lognormal, conditioned on [lo, hi] by redrawing the samples that fall outside
    mu, sigma = math.log(dist["median"]), float(dist["sigma"])
    out = np.empty(0)
    for _ in range(1000):
        x = np.exp(rng.normal(mu, sigma, max(2 * (n - len(out)), 16)))
        out = np.concatenate([out, x[(x >= lo) & (x <= hi)]])
        if len(out) >= n:
            return out[:n]
    raise ValueError(f"lognormal(median={dist['median']}, sigma={sigma}) almost never falls in "
                     f"[{lo}, {hi}]; widen the bounds")
