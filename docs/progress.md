# Progress

Status of the Scrap Shells project (6.8410 Shape Analysis). Last updated 2026-11-24.

For how the baseline works, see [baseline.md](baseline.md). For the stock of offcuts, see
[stock.md](stock.md). For what comes next, see [plan.md](plan.md).

## At a glance

| Milestone | Status | Where |
| --- | --- | --- |
| Proposal | Done (2026-10-13) | Proposal report and slide deck (not in the repo) |
| Baseline: PH-CPF reproduced on macOS | Done (2026-11-11) | [README.md](../README.md), [baseline.md](baseline.md) |
| Stage-by-stage web viewer | Done (2026-11-11) | [zincr.github.io/ScrapShell/web](https://zincr.github.io/ScrapShell/web/) |
| Checkpoint report | Done (Nov 2026) | Not in the repo |
| Extension: stock generator | Done (2026-11-22) | [stock.md](stock.md), [stock/](../stock/) |
| Extension: stock-driven panel layout | In progress: footprint model done (2026-11-24); plan revised the same day | [plan.md](plan.md), [sizing/](../sizing/) |
| Final report and presentation | Not started | |

## What works

The full PH-CPF pipeline runs natively on an Apple Silicon Mac. The official MATLAB code handles
curvature, the field optimization and integer rounding, with the official parameters. Our own code
replaces the Windows-only executables: cutting, the Poisson system and meshing in C++ on Directional,
and dualization and planarization in MATLAB.

**6 of the 8 models shipped with the official code reach the paper's 1% planarity target.** Every
panel ends within 1% planarity (paper Eq. 29), and the panel mesh stays within 0.72% of the bounding-box
diagonal (bbd) of the input surface.

| Model | Faces in | Singularities | Panels (hex + other) | Planarity max, before → after | Planarity mean | Distance max / bbd | Pipeline time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Metro | 1418 | 21 | 422 + 14 | 4.60% → 0.87% | 0.122% | 0.0072 | 19.7 s |
| Intersection | 1423 | 22 | 502 + 20 | 8.88% → 0.95% | 0.120% | 0.0050 | 19.6 s |
| Six | 1365 | 5 | 412 + 4 | 5.60% → 0.77% | 0.171% | 0.0055 | 11.3 s |
| Vase | 2404 | 0 | 394 + 0 | 2.14% → 0.85% | 0.251% | 0.0017 | 7.0 s |
| Sphere with cones | 1972 | 12 | 702 + 12 | 2.72% → 0.86% | 0.252% | 0.0029 | 11.4 s |
| C | 996 | 14 | 877 + 12 | 5.75% → 0.85% | 0.180% | 0.0062 | 12.4 s |

"Other" panels are the pentagons and heptagons around singularities. Pipeline time is the sum of the
stage timings on an Apple M5 Pro and excludes MATLAB startup. Planarization takes 55–66% of it on
Metro, Intersection, Six and C.

All numbers fall inside the ranges the paper reports (Fig. 24/25): 255–1449 panels, planarity max
0.3–1.0%, mean 0.03–0.49%, distance 0.001–0.015 × bbd. Our maxima sit just under 1% because the
planarizer stops as soon as the target is met.

To regenerate: `scripts/run_all.sh`, then `python scripts/results_table.py`. The per-run numbers are
in `web/data/<model>.meta.json`.

## What does not work

**Car and Three holes stop at integer rounding (stage 5).** Both are closed models with many
singularities: 52 and 32 on the last attempt. The official `BarrierPoissonSolve.m` cannot find an
integer solution without flipped triangles, even after the paper's continuation: the continuity
weight doubled 10 times (the official `runme.m` stops after 3). The viewer shows these models up to
the field stage. We plan to leave them as a documented limitation.

## Fixes made along the way

These are the changes we needed beyond replacing the executables. Details are in
[baseline.md](baseline.md#how-it-differs-from-the-official-release).

1. **Directional mesher bug.** Coincident isoline edges were merged only when they had the same
   orientation. An isoline running exactly along a mesh edge then broke the mesh (empty output). We
   patch a build-local copy of the header; the submodule is untouched.
2. **Planarizer normal collapse.** With the panel normal entering the planarity term directly, the
   growing planarity weight made shrinking a hard panel's normal to zero cheaper than flattening the
   panel. A few panels stalled at 5–7%. The planarity term now uses only the normal's direction.
3. **Robustness.** Each mesher call is capped at 120 s. Rounding retries go up to 10 continuation
   levels. `OPT.m` drops an option MATLAB R2026a rejects. The Windows-only `fast_sparse` MEX file is
   replaced by `sparse`.

## Known loose ends

- The comment at the top of [scripts/run_all.sh](../scripts/run_all.sh) says failed models are retried
  with `iter_opt` 1..3; [matlab/run_with_retries.m](../matlab/run_with_retries.m) actually goes to 10.
- `results/smoke/` holds debugging runs from the port. It is git-ignored and can be deleted.

## Extension so far

**The stock exists** ([stock.md](stock.md)). One deposit of 4000 scraps in meters,
[data/stock/deposit_s0.json](../data/stock/deposit_s0.json), shared by every model: 1000 each of
four shape families (remnant, shard, blob, board), all with the same size distribution (√area
lognormal, median 0.5 m). Every model is scaled to a 15 m span, so δ is 75 mm on every model. The
deposit holds 4–16× each model's surface area, so piece count never limits a layout.

**At baseline panel sizes, fit is the limit.** A quick fit preview (`python -m stock.plot`) finds
that only 7–38% of scraps can hold a typical baseline panel when it is centred on the scrap's
incenter; the exact test will find somewhat more. Metro is tightest: about 300 suitable scraps for
436 panels. Blobs and shards fit hexagonal panels best, boards worst. We keep the 0.5 m median so
that the stock genuinely shapes the layout.

**The footprint model exists** ([sizing/footprint.py](../sizing/footprint.py)). It recovers each
face's curvatures from the baseline's sizing ellipses and predicts the panel PH-CPF makes there: a
regular hexagon in the sizing metric with circumradius 0.70 matches the panels of all six finished
models (0.68–0.71 each) and their counts within 6%.

**It changed the plan.** The maximum element area η, not the curvature, sets the panel size on
65–100% of each finished model's surface. Scaling the metric by a factor σ ≤ 1, as first planned,
could therefore never let panels grow toward large scraps without raising η everywhere, and that
stretches panels (up to 5.7:1 at the 90th percentile), which the deposit fits poorly. The size field now sets η per face
instead ([plan.md](plan.md#why-the-cap-varies-per-face)).

**Not started yet:** the η hook in CPF, scrap capacity and the coupling, the exact fit test,
panel-to-scrap assignment and the feedback loop. [plan.md](plan.md) lays out the work.

## Log

- **2026-10-13** Proposal submitted: report and slide deck.
- **2026-11-11** PH-CPF baseline committed (`e9a17c9`), with the web viewer. The viewer is served
  from GitHub Pages (`26fcc55`).
- **Nov 2026** Checkpoint report submitted: baseline on 6 of 8 models.
- **2026-11-21** Extension plan drafted ([plan.md](plan.md)); this `docs/` folder added.
- **2026-11-21** Stock generator, steps 1–3: shape descriptors (`65a40ba`), four shape families
  (`a7e955f`), and one deposit in meters with every model at a 15 m span (`1ff7913`). Design
  changes along the way: one deposit for all models instead of one per model, equal family
  shares, and absolute sizes instead of sizes relative to each model.
- **2026-11-22** Stock figures and fit preview (`560dab6`): 7–38% of scraps hold a typical
  baseline panel. Decision: keep the 0.5 m median size. [stock.md](stock.md) added.
- **2026-11-24** Size field, step 1: the footprint model (`sizing/footprint.py`), with one panel
  constant fitted on the baseline. Decision: a per-face η replaces the size factor σ, because η caps
  the panel size on most of every surface and raising it for all faces stretches panels;
  [plan.md](plan.md) revised.
