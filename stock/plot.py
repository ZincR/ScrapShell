"""Figures of the deposit against the models' baseline panels, and a quick fit preview.

    python -m stock.plot                              # results/stock/*.png, fit_preview.json, summary
    python -m stock.plot --out docs/figures/stock     # the figures in docs/stock.md

Figures (every model at the common span of stock/place.py, so everything is in meters):
  deposit_sample.png  scraps of each family across the size range, and each model's median
                      baseline panel, to scale
  deposit_pile.png    the whole deposit laid out side by side, to scale
  sizes.png           scrap areas against each model's baseline panel areas
  fit_preview.png     how much of the deposit can hold a model's baseline panels

The fit preview is a quick estimate, not the fit test of docs/plan.md (Phase 2, item 5). For a
sample of each model's baseline panels it counts the scraps that pass three tests, each implying
the one after it:
  circle:   the scrap's inscribed circle is at least the panel's enclosing circle (coarse);
  incenter: the panel, with its enclosing circle centred on the scrap's incenter, lies inside the
            scrap at one of `rotations` angles (the exact test will also slide the panel, and may
            flip it, so this is a lower bound);
  area:     the scrap's area is at least the panel's (necessary, so an upper bound).
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import shapely  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from stock.generate import DEFAULT_DEPOSIT, load_deposit  # noqa: E402
from stock.model_ref import ROOT, load_model_ref  # noqa: E402
from stock.place import DEFAULT_SPAN, model_ids, placement  # noqa: E402
from stock.shapes import FAMILIES  # noqa: E402

OUT_DIR = ROOT / "results" / "stock"

# reference palette (dataviz skill): chart surface, ink, and two categorical slots
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
SCRAP, PANEL = "#2a78d6", "#eb6834"          # slot 1 blue: scraps; slot 2 orange: baseline panels


def _style():
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 10, "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": AXIS,
        "axes.linewidth": 0.8, "xtick.color": AXIS, "ytick.color": AXIS,
        "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
        "axes.spines.top": False, "axes.spines.right": False, "grid.color": GRID,
        "grid.linewidth": 0.8, "legend.frameon": False, "axes.titlelocation": "left",
        "axes.titlesize": 12, "axes.titleweight": "bold",
    })


def _title(fig, title: str, subtitle: str):
    fig.text(0.012, 0.985, title, ha="left", va="top", fontsize=13, fontweight="bold", color=INK)
    fig.text(0.012, 0.945, subtitle, ha="left", va="top", fontsize=9.5, color=INK2)


def _scale_bar(ax, x, y, length: float, label: str):
    ax.plot([x, x + length], [y, y], color=INK2, lw=1.5, solid_capstyle="butt")
    for xx in (x, x + length):
        ax.plot([xx, xx], [y - 0.06 * length, y + 0.06 * length], color=INK2, lw=1)
    ax.text(x + length / 2, y - 0.12 * length, label, ha="center", va="top", fontsize=9, color=INK2)


def _completed(models):
    """(model id, ModelRef, Placement) for the models whose baseline produced panels."""
    out = []
    for m in models:
        ref = load_model_ref(m)
        if ref.panels:
            out.append((m, ref, placement(m)))
    return out


# --- figures ----------------------------------------------------------------------------------

def deposit_sample(deposit, models, path: Path, per_family: int = 14):
    """Scraps of each family at evenly spaced size ranks, and each model's median panel, to scale."""
    gap, row_gap = 0.12, 0.45
    rows = []
    for fam in FAMILIES:
        scraps = sorted((s for s in deposit["scraps"] if s["family"] == fam), key=lambda s: s["size"])
        pick = np.unique(np.linspace(0, len(scraps) - 1, per_family).round().astype(int))
        rows.append((fam, [np.array(scraps[i]["polygon"]) for i in pick], SCRAP))
    panels, names = [], []
    for m, ref, pl in _completed(models):
        k = int(np.argmin(np.abs(ref.panel_areas - np.median(ref.panel_areas))))
        panels.append(ref.panels[k] * pl.meters_per_unit)
        names.append(ref.label)
    rows.append(("baseline panel", panels, PANEL))

    fig, ax = plt.subplots(figsize=(13, 7.2))
    y, width = 0.0, 0.0
    for label, polys, color in rows:
        h = max(p[:, 1].max() - p[:, 1].min() for p in polys)
        y -= h
        x, verts = 0.0, []
        for j, p in enumerate(polys):
            w = p[:, 0].max() - p[:, 0].min()
            slot = max(w, 1.15) if label == "baseline panel" else w     # room for the model's name
            q = p - [p[:, 0].min(), p[:, 1].min()] + [x + (slot - w) / 2, y]
            verts.append(q)
            if label == "baseline panel":
                ax.text(x + slot / 2, y - 0.08, names[j], ha="center", va="top", fontsize=8.5,
                        color=INK2)
            x += slot + gap
        width = max(width, x)
        ax.add_collection(PolyCollection(verts, facecolors=color, edgecolors=SURFACE, linewidths=0.8))
        ax.text(-0.25, y + h / 2, label, ha="right", va="center", fontsize=10, color=INK)
        y -= row_gap
    _scale_bar(ax, width - 1.0, y + 0.1, 1.0, "1 m")
    ax.set_xlim(-2.2, width + 0.1)
    ax.set_ylim(y - 0.2, 0.15)
    ax.set_aspect("equal")
    ax.axis("off")
    _title(fig, "The deposit, to scale",
           f"{per_family} scraps per family from smallest to largest; bottom row: each model's median "
           f"baseline panel at a {DEFAULT_SPAN:g} m span")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.02)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def deposit_pile(deposit, path: Path):
    """All scraps packed in rows (tallest first), to scale."""
    polys = [np.array(s["polygon"]) for s in deposit["scraps"]]
    boxes = [(p[:, 0].max() - p[:, 0].min(), p[:, 1].max() - p[:, 1].min()) for p in polys]
    total = deposit["summary"]["totalArea"]
    row_w, gap = math.sqrt(total / 0.6 * 1.45), 0.06     # packing density ~0.6, layout 1.45:1
    order = sorted(range(len(polys)), key=lambda i: -boxes[i][1])
    verts, x, y, row_h = [], 0.0, 0.0, 0.0
    for i in order:
        w, h = boxes[i]
        if x + w > row_w:
            x, y, row_h = 0.0, y - row_h - gap, 0.0
        row_h = max(row_h, h)
        p = polys[i]
        verts.append(p - [p[:, 0].min(), p[:, 1].max()] + [x, y])
        x += w + gap
    fig, ax = plt.subplots(figsize=(11, 7.6))
    ax.add_collection(PolyCollection(verts, facecolors=SCRAP, edgecolors=SURFACE, linewidths=0.3))
    bottom = y - row_h
    _scale_bar(ax, row_w - 5.0, bottom - 1.0, 5.0, "5 m")
    ax.set_xlim(-0.5, row_w + 0.5)
    ax.set_ylim(bottom - 2.8, 0.5)
    ax.set_aspect("equal")
    ax.axis("off")
    s = deposit["summary"]
    _title(fig, "The whole deposit",
           f"{s['count']} scraps, {s['totalArea']:,.0f} m², 1000 per family; packed side by side, "
           f"tallest first")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.91, bottom=0.01)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def sizes_chart(deposit, models, path: Path):
    """Range plot: scrap areas against each model's baseline panel areas (m^2, log scale)."""
    rows = [("Deposit scraps", np.array([s["area"] for s in deposit["scraps"]]), SCRAP)]
    for m, ref, pl in _completed(models):
        rows.append((ref.label, ref.panel_areas * pl.meters_per_unit ** 2, PANEL))
    fig, ax = plt.subplots(figsize=(9, 0.55 * len(rows) + 1.9))
    for k, (label, a, color) in enumerate(rows):
        yk = -k
        p10, p50, p90 = np.percentile(a, [10, 50, 90])
        ax.plot([a.min(), a.max()], [yk, yk], color=AXIS, lw=1, zorder=1)
        ax.plot([p10, p90], [yk, yk], color=color, lw=5, solid_capstyle="round", zorder=2)
        ax.plot([p50], [yk], "o", ms=9, color=color, mec=SURFACE, mew=2, zorder=3)
    ax.set_yticks([-k for k in range(len(rows))], [r[0] for r in rows])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xscale("log")
    ticks = [0.02, 0.05, 0.1, 0.2, 0.5, 1, 2]
    ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    ax.set_xlabel("Area (m²)")
    ax.set_ylim(-len(rows) + 0.4, 0.6)
    ax.legend(handles=[Line2D([], [], color=SCRAP, lw=5, label="Scraps"),
                       Line2D([], [], color=PANEL, lw=5, label="Baseline panels")],
              loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, fontsize=9)
    _title(fig, "Scrap areas against baseline panel areas",
           f"Bar: 10th to 90th percentile; dot: median; line: full range. Models at a {DEFAULT_SPAN:g} m span.")
    fig.subplots_adjust(left=0.2, right=0.97, top=0.78, bottom=0.2)
    fig.savefig(path, dpi=150)
    plt.close(fig)


# --- fit preview ------------------------------------------------------------------------------

def fit_preview(deposit, models, n_panels: int = 16, rotations: int = 36, seed: int = 0) -> list[dict]:
    """For each completed model, the median (over a sample of its baseline panels) fraction of the
    deposit passing the circle, incenter and area tests (see the module docstring)."""
    scraps = deposit["scraps"]
    polys = [shapely.Polygon(s["polygon"]) for s in scraps]
    shapely.prepare(polys)
    area = np.array([s["area"] for s in scraps])
    inradius = np.array([s["inradius"] for s in scraps])
    family = np.array([s["family"] for s in scraps])
    angles = np.linspace(0, 360, rotations, endpoint=False)
    rng = np.random.default_rng(seed)
    results = []
    for m, ref, pl in _completed(models):
        pick = rng.choice(len(ref.panels), size=min(n_panels, len(ref.panels)), replace=False)
        per_panel = []
        for i in pick:
            panel = shapely.Polygon(ref.panels[i] * pl.meters_per_unit)
            circle_center = shapely.minimum_bounding_circle(panel).centroid
            radius = shapely.minimum_bounding_radius(panel)
            panel = shapely.affinity.translate(panel, -circle_center.x, -circle_center.y)
            turned = [shapely.affinity.rotate(panel, a, origin=(0, 0)) for a in angles]
            circle_ok = inradius >= radius
            fits = circle_ok.copy()
            for j in np.flatnonzero((area >= panel.area) & ~circle_ok):
                fits[j] = shapely.contains(polys[j], turned).any()
            per_panel.append({"area": np.mean(area >= panel.area), "incenter": np.mean(fits),
                              "circle": np.mean(circle_ok),
                              **{f"family:{f}": np.mean(fits[family == f]) for f in FAMILIES}})
        med = {k: float(np.median([p[k] for p in per_panel])) for k in per_panel[0]}
        results.append({"model": m, "label": ref.label, "panels": len(pick),
                        "area": med["area"], "incenter": med["incenter"], "circle": med["circle"],
                        "byFamily": {f: med[f"family:{f}"] for f in FAMILIES}})
    return results


def fit_preview_chart(results, path: Path, rotations: int = 36):
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12, 0.5 * len(results) + 2.6),
                                 gridspec_kw={"width_ratios": [3, 2], "wspace": 0.45})
    for k, r in enumerate(results):
        lo, mid, hi = 100 * r["circle"], 100 * r["incenter"], 100 * r["area"]
        ax.plot([lo, hi], [-k, -k], color=AXIS, lw=1.5, zorder=1)
        ax.plot([lo], [-k], marker="<", ms=7, color=MUTED, ls="none", zorder=2)
        ax.plot([hi], [-k], marker=">", ms=7, color=MUTED, ls="none", zorder=2)
        ax.plot([mid], [-k], "o", ms=9, color=SCRAP, mec=SURFACE, mew=2, zorder=3)
        ax.text(mid, -k + 0.32, f"{mid:.0f}%", ha="center", va="bottom", fontsize=8.5, color=INK2)
    ax.set_yticks([-k for k in range(len(results))], [r["label"] for r in results])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 100)
    ax.set_ylim(-len(results) + 0.4, 0.8)
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    ax.set_xlabel("Share of the deposit that can hold a typical panel (%)")
    ax.legend(handles=[
        Line2D([], [], marker="o", ms=8, color=SCRAP, ls="none", label=f"Fits, centred on the incenter ({rotations} rotations)"),
        Line2D([], [], marker=">", ms=7, color=MUTED, ls="none", label="Upper bound: scrap area ≥ panel area"),
        Line2D([], [], marker="<", ms=7, color=MUTED, ls="none", label="Lower bound: inscribed circle holds the panel"),
    ], loc="upper left", bbox_to_anchor=(0.0, -0.2 - 0.02 * len(results)), fontsize=8.5, ncol=1)

    fams = list(FAMILIES)
    share = [100 * np.mean([r["byFamily"][f] for r in results]) for f in fams]
    bx.barh(range(len(fams))[::-1], share, height=0.5, color=SCRAP)
    for k, v in zip(range(len(fams))[::-1], share):
        bx.text(v + 1.5, k, f"{v:.0f}%", va="center", fontsize=8.5, color=INK2)
    bx.set_yticks(range(len(fams))[::-1], fams)
    bx.tick_params(axis="y", length=0)
    bx.spines["left"].set_visible(False)
    bx.set_xlim(0, 100)
    bx.grid(axis="x")
    bx.set_axisbelow(True)
    bx.set_xlabel("Fits, centred on the incenter (%), mean over models")
    bx.set_title("By family", fontsize=10.5, color=INK, pad=8)
    _title(fig, "Fit preview: how much of the deposit can hold a model's baseline panels",
           f"Median over {results[0]['panels']} sampled panels per model, at a {DEFAULT_SPAN:g} m span. "
           f"A quick estimate; the exact fit test also slides and may flip panels.")
    fig.subplots_adjust(left=0.14, right=0.97, top=0.8, bottom=0.36)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m stock.plot",
                                 description="Figures of the deposit against the baseline panels.")
    ap.add_argument("--deposit", type=Path, default=DEFAULT_DEPOSIT, help="deposit file")
    ap.add_argument("--out", type=Path, default=OUT_DIR, help="output folder (default results/stock)")
    ap.add_argument("--models", nargs="+", help="models to include (default: all in the manifest)")
    ap.add_argument("--panels", type=int, default=16, help="baseline panels sampled per model for the fit preview")
    ap.add_argument("--rotations", type=int, default=36, help="rotations tried in the fit preview")
    args = ap.parse_args(argv)

    _style()
    deposit = load_deposit(args.deposit)
    models = args.models or model_ids()
    args.out.mkdir(parents=True, exist_ok=True)
    deposit_sample(deposit, models, args.out / "deposit_sample.png")
    deposit_pile(deposit, args.out / "deposit_pile.png")
    sizes_chart(deposit, models, args.out / "sizes.png")
    results = fit_preview(deposit, models, args.panels, args.rotations)
    fit_preview_chart(results, args.out / "fit_preview.png", args.rotations)
    (args.out / "fit_preview.json").write_text(json.dumps(results, indent=2) + "\n")

    print(f"figures in {args.out}")
    print(f"{'model':16s} {'circle':>7s} {'incenter':>9s} {'area':>6s}   by family (incenter)")
    for r in results:
        fam = "  ".join(f"{f} {100 * v:3.0f}%" for f, v in r["byFamily"].items())
        print(f"{r['label']:16s} {100 * r['circle']:6.0f}% {100 * r['incenter']:8.0f}% "
              f"{100 * r['area']:5.0f}%   {fam}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
