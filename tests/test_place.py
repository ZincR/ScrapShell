"""Tests for stock/place.py: models at a common span, and the deposit in model units."""
import numpy as np
import pytest
import shapely

from stock.generate import dumps_deposit, generate_deposit
from stock.model_ref import load_model_ref
from stock.place import DEFAULT_SPAN, main, model_ids, place, placement


def test_metro_at_the_default_span():
    pl = placement("metro")
    ref = load_model_ref("metro")
    assert pl.span == DEFAULT_SPAN == 15.0
    assert pl.meters_per_unit == pytest.approx(15.0 / ref.bbd)
    assert pl.surface_area == pytest.approx(ref.surface_area * (15.0 / ref.bbd) ** 2)
    assert pl.delta == pytest.approx(0.075)                       # 0.5% of 15 m
    assert pl.n_panels == 436
    assert pl.panel_area_median == pytest.approx(ref.a_ref * pl.meters_per_unit ** 2)


@pytest.mark.parametrize("model", model_ids())
def test_every_model_gets_the_same_span(model):
    """Including the models the baseline did not finish: placement needs only the surface."""
    pl = placement(model, span=20.0)
    assert load_model_ref(model).bbd * pl.meters_per_unit == pytest.approx(20.0)
    assert pl.delta == pytest.approx(0.1)
    assert (pl.panel_area_median is None) == (pl.n_panels == 0)


def test_place_converts_meters_to_model_units():
    deposit = generate_deposit({"n_scraps": 8}, seed=0)
    pl = placement("vase")
    outlines = place(deposit, pl)
    assert len(outlines) == 8
    for s, q in zip(deposit["scraps"], outlines):
        assert np.allclose(q * pl.meters_per_unit, s["polygon"])
        assert shapely.Polygon(q).area == pytest.approx(s["area"] / pl.meters_per_unit ** 2, rel=1e-6)


def test_bad_span():
    with pytest.raises(ValueError, match="positive"):
        placement("metro", span=0)


def test_command_line(tmp_path, capsys):
    deposit = tmp_path / "d.json"
    deposit.write_text(dumps_deposit(generate_deposit({"n_scraps": 8}, seed=0)))
    assert main(["--deposit", str(deposit)]) == 0
    out = capsys.readouterr().out
    assert all(m in out for m in model_ids())
    assert "8 scraps" in out
