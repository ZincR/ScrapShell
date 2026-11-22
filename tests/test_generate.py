"""Tests for stock/generate.py: the deposit document and the command line."""
import json
import math
from collections import Counter

import numpy as np
import pytest
from shapely.geometry import Point, Polygon

from stock.describe import describe
from stock.generate import (DEFAULTS, dumps_deposit, generate_deposit, load_deposit, main,
                            resolve_config, sig)
from stock.shapes import FAMILIES, check

DESCRIPTOR_KEYS = ["area", "perimeter", "nVertices", "inradius", "incenter", "circumradius",
                   "length", "width", "elongation", "axisAngle", "convexity"]
N = 200


@pytest.fixture(scope="module")
def doc():
    return generate_deposit({"n_scraps": N}, seed=0)


def test_header(doc):
    assert doc["format"] == "scrapshell-stock/2"
    assert doc["name"] == "deposit_s0" and doc["units"] == "m" and doc["seed"] == 0
    assert doc["generator"]["n_scraps"] == N
    assert doc["generator"]["sizes"] == DEFAULTS["sizes"]
    assert doc["generator"]["families"] == list(FAMILIES)
    assert doc["summary"]["count"] == len(doc["scraps"]) == N


def test_equal_family_shares(doc):
    assert Counter(s["family"] for s in doc["scraps"]) == {f: N // 4 for f in FAMILIES}


def test_every_family_gets_the_same_sizes(doc):
    """Sorted by size, each consecutive group of four holds one scrap of each family."""
    by_size = sorted(doc["scraps"], key=lambda s: s["size"])
    for k in range(0, N, 4):
        assert {s["family"] for s in by_size[k:k + 4]} == set(FAMILIES)


def test_ids_are_shuffled(doc):
    """Neither family nor size runs in blocks along the ids."""
    first = doc["scraps"][:40]
    assert {s["family"] for s in first} == set(FAMILIES)
    ids, size = np.arange(N), np.array([s["size"] for s in doc["scraps"]])
    assert abs(np.corrcoef(ids, size)[0, 1]) < 0.3


def test_sizes(doc):
    lo, hi = DEFAULTS["sizes"]["min"], DEFAULTS["sizes"]["max"]
    for s in doc["scraps"]:
        assert lo * (1 - 1e-6) <= s["size"] <= hi * (1 + 1e-6)
        assert s["size"] == pytest.approx(math.sqrt(s["area"]), rel=1e-6)
    assert np.median([s["size"] for s in doc["scraps"]]) == pytest.approx(0.5, rel=0.1)


def test_scraps(doc):
    for i, s in enumerate(doc["scraps"]):
        assert s["id"] == i
        poly = np.array(s["polygon"])
        assert check(poly) is None
        assert {k: s[k] for k in DESCRIPTOR_KEYS} == sig(describe(poly).to_dict())
        p = Polygon(poly)                    # canonical frame, in meters
        assert p.contains(Point(0, 0))
        assert p.exterior.distance(Point(0, 0)) == pytest.approx(s["inradius"], rel=1e-3)


def test_reproducible_and_seed_dependent():
    a = dumps_deposit(generate_deposit({"n_scraps": 40}, seed=3))
    assert a == dumps_deposit(generate_deposit({"n_scraps": 40}, seed=3))
    assert a != dumps_deposit(generate_deposit({"n_scraps": 40}, seed=4))


def test_family_parameters_only_change_that_family():
    a = generate_deposit({"n_scraps": 80}, seed=5)
    b = generate_deposit({"n_scraps": 80, "params": {"blob": {"spikiness": 0.1}}}, seed=5)
    for sa, sb in zip(a["scraps"], b["scraps"]):
        assert sa["family"] == sb["family"] and sa["size"] == pytest.approx(sb["size"], rel=1e-5)
        if sa["family"] != "blob":
            assert sa["polygon"] == sb["polygon"]
    assert any(sa["polygon"] != sb["polygon"] for sa, sb in zip(a["scraps"], b["scraps"])
               if sa["family"] == "blob")


@pytest.mark.parametrize("config, message", [
    ({"n_scraps": 10}, "multiple of 4"),
    ({"n_scraps": 0}, "multiple of 4"),
    ({"n_scraps": 4.0}, "multiple of 4"),
    ({"sizes": {"kind": "gamma"}}, "unknown size distribution"),
    ({"params": {"plank": {}}}, "unknown families"),
    ({"colour": "oak"}, "unknown config keys"),
])
def test_bad_config(config, message):
    with pytest.raises(ValueError, match=message):
        resolve_config(config)


def test_limits_are_filled_in():
    assert resolve_config({"limits": {"min_angle": 30.0}})["limits"]["min_angle"] == 30.0
    assert resolve_config()["limits"]["max_vertices"] == 16


def test_dumps_round_trip(doc, tmp_path):
    text = dumps_deposit(doc)
    assert json.loads(text) == doc
    assert sum(line.startswith('    {"id":') for line in text.splitlines()) == N
    path = tmp_path / "d.json"
    path.write_text(text)
    assert load_deposit(path) == doc


def test_load_rejects_other_formats(tmp_path):
    path = tmp_path / "old.json"
    path.write_text(json.dumps({"format": "scrapshell-stock/1"}))
    with pytest.raises(ValueError, match="not a scrapshell-stock/2"):
        load_deposit(path)


def test_command_line(tmp_path, capsys):
    out, cfg = tmp_path / "deposit.json", tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"sizes": {"kind": "uniform", "min": 0.2, "max": 0.4}}))
    assert main(["--seed", "2", "--n-scraps", "20", "--config", str(cfg), "-o", str(out)]) == 0
    expected = generate_deposit({"n_scraps": 20, "sizes": {"kind": "uniform", "min": 0.2, "max": 0.4}}, 2)
    assert out.read_text() == dumps_deposit(expected)
    assert "20 scraps" in capsys.readouterr().out


def test_command_line_error(tmp_path):
    with pytest.raises(SystemExit) as exc:
        main(["--n-scraps", "10", "-o", str(tmp_path / "x.json")])
    assert exc.value.code == 2
    assert not (tmp_path / "x.json").exists()
