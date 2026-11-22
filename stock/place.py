"""Put the deposit on a model: give the model a physical size, and convert meters to model units.

The deposit is in meters (stock/generate.py); the model meshes are in arbitrary units. Every model
gets the same span: its bounding-box diagonal is scaled to `span` meters (default 15 m). PH-CPF's
tolerances are fractions of the model's size (delta = 0.5% of the diagonal, as in runme.m), so the
baseline's results simply scale with the model, and with a common span delta is the same distance
in meters on every model (0.5% of 15 m = 75 mm).

    python -m stock.place                 # each model's size at the default span, and the stock's coverage
    python -m stock.place --span 20
"""
from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from stock.generate import DEFAULT_DEPOSIT, load_deposit
from stock.model_ref import WEB_DATA, load_model_ref

DEFAULT_SPAN = 15.0      # m: the bounding-box diagonal given to every model
DEFAULT_DELTA = 0.005    # PH-CPF delta: maximum panel-to-surface distance, fraction of the diagonal


@dataclass(frozen=True)
class Placement:
    """A model at its physical size."""
    model: str
    span: float                 # m: bounding-box diagonal
    meters_per_unit: float      # length of one model unit, in m
    surface_area: float         # m^2
    delta: float                # m: PH-CPF's maximum panel-to-surface distance at this span
    n_panels: int               # baseline panels (0 if the baseline did not planarize this model)
    panel_area_median: float | None   # m^2: the baseline's median panel

    def to_model_units(self, coords_m) -> np.ndarray:
        """Coordinates in meters -> the model's units."""
        return np.asarray(coords_m, dtype=float) / self.meters_per_unit

    def to_dict(self) -> dict:
        return {"model": self.model, "span": self.span, "metersPerUnit": self.meters_per_unit,
                "surfaceArea": self.surface_area, "delta": self.delta, "nPanels": self.n_panels,
                "panelAreaMedian": self.panel_area_median}


def placement(model: str, span: float = DEFAULT_SPAN, delta: float = DEFAULT_DELTA,
              data_dir: Path = WEB_DATA) -> Placement:
    """The model scaled so that its bounding-box diagonal is `span` meters."""
    if span <= 0:
        raise ValueError("span must be positive")
    ref = load_model_ref(model, data_dir)
    s = span / ref.bbd
    return Placement(model=model, span=span, meters_per_unit=s, surface_area=ref.surface_area * s**2,
                     delta=delta * span, n_panels=len(ref.panels),
                     panel_area_median=ref.a_ref * s**2 if ref.panels else None)


def place(deposit: Mapping, pl: Placement) -> list[np.ndarray]:
    """The deposit's scrap outlines in the model's units, in the order of deposit["scraps"]."""
    return [pl.to_model_units(s["polygon"]) for s in deposit["scraps"]]


def model_ids(data_dir: Path = WEB_DATA) -> list[str]:
    """Models in the viewer's manifest, in its order."""
    return [m["id"] for m in json.loads((Path(data_dir) / "manifest.json").read_text())["models"]]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m stock.place",
                                 description="Each model at a common span, against the deposit.")
    ap.add_argument("--span", type=float, default=DEFAULT_SPAN, help=f"bounding-box diagonal in m (default {DEFAULT_SPAN:g})")
    ap.add_argument("--deposit", type=Path, default=DEFAULT_DEPOSIT, help="deposit file (default data/stock/deposit_s0.json)")
    args = ap.parse_args(argv)

    deposit = load_deposit(args.deposit) if args.deposit.exists() else None
    if deposit:
        s = deposit["summary"]
        print(f"deposit {args.deposit.name}: {s['count']} scraps, {s['totalArea']:g} m^2, "
              f"size median {s['size']['median']:g} m")
    else:
        print(f"no deposit at {args.deposit}; run python -m stock.generate first")
    print(f"{'model':14s} {'m/unit':>8s} {'area m^2':>9s} {'delta mm':>9s} {'panels':>7s} "
          f"{'panel m^2':>10s} {'stock/area':>11s} {'scraps/panel':>13s}")
    for m in model_ids():
        pl = placement(m, args.span)
        panel = f"{pl.panel_area_median:.3f}" if pl.panel_area_median else "-"
        cover = f"{deposit['summary']['totalArea'] / pl.surface_area:.1f}x" if deposit else "-"
        per = f"{deposit['summary']['count'] / pl.n_panels:.1f}" if deposit and pl.n_panels else "-"
        print(f"{m:14s} {pl.meters_per_unit:8.4g} {pl.surface_area:9.1f} {1000 * pl.delta:9.0f} "
              f"{pl.n_panels or '-':>7} {panel:>10s} {cover:>11s} {per:>13s}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
