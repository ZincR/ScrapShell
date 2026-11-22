"""Tests for stock/plot.py: the figures render, and the fit preview's bounds are ordered."""
import json

from stock.generate import dumps_deposit, generate_deposit
from stock.plot import main
from stock.shapes import FAMILIES

FIGURES = ["deposit_sample.png", "deposit_pile.png", "sizes.png", "fit_preview.png"]


def test_figures_and_fit_preview(tmp_path, capsys):
    deposit = tmp_path / "deposit.json"
    deposit.write_text(dumps_deposit(generate_deposit({"n_scraps": 40}, seed=0)))
    out = tmp_path / "figs"
    assert main(["--deposit", str(deposit), "--out", str(out), "--models", "metro", "intersection",
                 "car_r", "--panels", "3", "--rotations", "12"]) == 0
    for name in FIGURES:
        assert (out / name).stat().st_size > 10_000, name
    results = json.loads((out / "fit_preview.json").read_text())
    assert [r["model"] for r in results] == ["metro", "intersection"]     # car_r has no panels
    for r in results:
        # each test implies the next: circle -> incenter -> area
        assert 0 <= r["circle"] <= r["incenter"] <= r["area"] <= 1
        assert set(r["byFamily"]) == set(FAMILIES)
        assert all(0 <= v <= 1 for v in r["byFamily"].values())
    assert "Metro" in capsys.readouterr().out
