"""Assign panels to scraps, one panel per scrap (docs/plan.md, Phase 2, item 6).

    python -m matching.assign metro            # baseline panels against the deposit
    python -m matching.assign metro_eta_4x_smooth

The cost of a pair is the waste, scrap area minus panel area; a pair that does not fit costs more
than any assignment of fitting pairs, so `linear_sum_assignment` first fits as many panels as
possible and then wastes as little as it can. Outputs the fit rate (share of panels assigned to a
scrap they fit inside), the scraps used, utilization (panel area / area of the scraps used), the
panels that fit nowhere, and each panel's placement in its scrap; written to
results/matching/<model>.json.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from shapely.geometry import Polygon

from matching.fit import ROTATIONS, FitGraph, fit_graph, load_deposit, model_panels
from stock.model_ref import ROOT

OUT_DIR = ROOT / "results" / "matching"


@dataclass
class Assignment:
    scrap: np.ndarray           # (n_panels,) scrap index, -1 where the panel fits nowhere
    panel_area: np.ndarray      # (n_panels,) m^2
    scrap_area: np.ndarray      # (n_scraps,) m^2

    @property
    def n_panels(self) -> int:
        return len(self.scrap)

    @property
    def assigned(self) -> np.ndarray:
        return self.scrap >= 0

    @property
    def fit_rate(self) -> float:
        return float(self.assigned.mean())

    @property
    def utilization(self) -> float:
        used = self.scrap[self.assigned]
        return float(self.panel_area[self.assigned].sum() / self.scrap_area[used].sum()) if len(used) else 0.0


def assign(fits: np.ndarray, panel_area: np.ndarray, scrap_area: np.ndarray) -> Assignment:
    """Most panels assigned to scraps they fit (fits: (n_panels, n_scraps) bool), then least waste."""
    cols = np.flatnonzero(fits.any(axis=0))
    scrap = np.full(len(panel_area), -1)
    if len(cols):
        waste = scrap_area[cols][None, :] - panel_area[:, None]
        big = float(len(panel_area) * scrap_area.max() + 1.0)
        cost = np.where(fits[:, cols], waste, big)
        r, c = linear_sum_assignment(cost)
        ok = fits[r, cols[c]]
        scrap[r[ok]] = cols[c[ok]]
    return Assignment(scrap=scrap, panel_area=np.asarray(panel_area, float), scrap_area=np.asarray(scrap_area, float))


def max_matching(fits: np.ndarray) -> int:
    """Size of a maximum matching (a check on assign)."""
    m = maximum_bipartite_matching(csr_matrix(fits.astype(np.int8)), perm_type="column")
    return int((m >= 0).sum())


def run(model: str, rotations: int = ROTATIONS, workers: int | None = None) -> dict:
    t0 = time.time()
    panels = model_panels(model)
    deposit = load_deposit()
    scraps = deposit["scraps"]
    scrap_area = np.array([s["area"] for s in scraps])
    panel_area = np.array([Polygon(P).area for P in panels])
    graph: FitGraph = fit_graph(panels, deposit, rotations, workers)
    t_fit = time.time() - t0
    fits = graph.fits()
    a = assign(fits, panel_area, scrap_area)
    assert a.assigned.sum() == max_matching(fits)

    # the inscribed-circle rule alone (sure fits only), for comparison
    circle = assign(_sure_fits(panels, scraps), panel_area, scrap_area)

    used = a.scrap[a.assigned]
    families = [scraps[j]["family"] for j in used]
    summary = {
        "model": model,
        "nPanels": a.n_panels,
        "nScraps": len(scraps),
        "rotations": rotations,
        "pairsSure": graph.n_sure,
        "pairsTested": graph.n_tested,
        "pairsTestedFit": graph.n_tested_fit,
        "scrapsFittingAny": int(fits.any(axis=0).sum()),
        "panelsFittingAny": int(fits.any(axis=1).sum()),
        "assigned": int(a.assigned.sum()),
        "fitRate": a.fit_rate,
        "utilization": a.utilization,
        "panelAreaMedian": float(np.median(panel_area)),
        "usedByFamily": {f: families.count(f) for f in sorted(set(s["family"] for s in scraps))},
        "inscribedCircleRule": {"assigned": int(circle.assigned.sum()), "fitRate": circle.fit_rate},
        "timeFit": t_fit,
    }
    placements = []
    for i in range(a.n_panels):
        j = int(a.scrap[i])
        if j >= 0:
            angle, t = graph.placements[i][j]
            placements.append({"panel": i, "scrap": j, "angle": float(angle), "t": [float(t[0]), float(t[1])]})
        else:
            placements.append({"panel": i, "scrap": -1})
    return {"summary": summary, "assignment": placements}


def _sure_fits(panels, scraps) -> np.ndarray:
    from matching.fit import describe_panel
    S_in = np.array([s["inradius"] for s in scraps])
    R = np.array([describe_panel(P).circumradius for P in panels])
    return S_in[None, :] >= R[:, None]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m matching.assign",
                                 description="Fit a model's panels into the deposit's scraps and assign them.")
    ap.add_argument("model", help="viewer data id, e.g. metro or metro_eta_4x_smooth")
    ap.add_argument("--rotations", type=int, default=ROTATIONS)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--out", type=Path, default=None, help="default results/matching/<model>.json")
    args = ap.parse_args(argv)
    res = run(args.model, args.rotations, args.workers)
    out = args.out or OUT_DIR / f"{args.model}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1))
    s = res["summary"]
    print(f"{s['model']}: {s['assigned']} of {s['nPanels']} panels fit a scrap ({s['fitRate']:.1%}); "
          f"utilization {s['utilization']:.1%}; {s['scrapsFittingAny']} scraps fit some panel, "
          f"{s['panelsFittingAny']} panels fit some scrap; inscribed-circle rule alone: "
          f"{s['inscribedCircleRule']['assigned']} ({s['inscribedCircleRule']['fitRate']:.1%}); "
          f"used by family {s['usedByFamily']}; fit test {s['timeFit']:.0f} s; {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
