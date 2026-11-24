# ScrapShell

6.8410 (Shape Analysis) course project: building discretized shells from a stock of irregular flat
offcuts. This repository currently holds the **baseline**: a working macOS (Apple Silicon)
reproduction of PH-CPF, *Planar Hexagonal Meshing using Coordinate Power Fields* (Pluta, Edelstein,
Vaxman, Ben-Chen, ACM TOG 40(4), 2021), plus a web viewer that replays every stage of the pipeline
on eight example shells.

## Docs

- [docs/progress.md](docs/progress.md): where the project stands, results so far, and a dated log.
- [docs/baseline.md](docs/baseline.md): how the PH-CPF baseline works, stage by stage, and where each
  piece lives.
- [docs/stock.md](docs/stock.md): the stock of offcuts: what is in it, how it is generated and
  placed on a model, its file format, and how well it fits the baseline's panels.
- [docs/plan.md](docs/plan.md): the plan for the stock-driven extension, experiments and timeline.

## Interactive viewer

**[Open the viewer](https://zincr.github.io/ScrapShell/web/)** (hosted on GitHub Pages from
[web/index.html](web/index.html)). It replays all eight stages on every model: curvature regions,
guiding field, the field optimization iteration by iteration, the parameterization, the triangle mesh,
the hexagonal dual and the planarization passes, with planarity colours and per-model numbers. A
"Sizing checks" group holds the extension's variant runs of Metro (per-face η; see
[docs/plan.md](docs/plan.md), Phase 1, item 4).

To run it locally, serve the folder (browsers refuse to load `web/data/` from a `file://` page):

```bash
scripts/serve.sh        # http://localhost:8765/   (or: python3 -m http.server -d web)
```

The empty `.nojekyll` file at the repo root makes GitHub Pages serve the files as they are, instead
of running Jekyll over the whole repository (which fails on Markdown inside the Directional submodule).

## Pipeline

The structure follows the official `CPF.solve` (`uv2mesh` -> `cut_and_comb` -> `param_and_mesh`
-> `planarize`); only the Windows-only executables are replaced.

```
official MATLAB              phcpf_mesher setup        official MATLAB          phcpf_mesher mesh     ScrapShell MATLAB
CPF field (LM)          -->  (MeshCutter +        --> BarrierPoissonSolve.m --> (Mesher)          --> dualize, planarize,
matching + singularities     Parametrization):        integer rounding with     isoline mesh of       metrics, export to
(compute_matching_and_sing)  cut, comb, Poisson       injectivity barrier       the rounded           web/data/*.json
                             system export                                      functions
```

| Stage | What | Where |
| --- | --- | --- |
| 1 | Input triangle mesh | `third_party/PH-CPF/data/*.off` |
| 2 | Curvature regions (elliptic / parabolic / hyperbolic / planar) | official `SHAPEOP.guiding_field2` |
| 3 | Guiding 2-RoSy field and sizing ellipses (metric g) | official `CPF.set_param` |
| 4 | Coordinate power field U, V; matching and singularities | official `CPF.opt_power_cfs`, `compute_matching_and_sing` |
| 5 | Seamless parameterization h1..h3 | `phcpf_mesher setup` + official `BarrierPoissonSolve.m` |
| 6 | Primal triangle-dominant mesh from integer isolines | `phcpf_mesher mesh` |
| 7 | Barycentric hexagonal dual | `matlab/dualize_bary.m` |
| 8 | Planarized hex mesh | `matlab/planarize_hex.m` |

## Setup (macOS, Apple Silicon)

1. **MATLAB** R2024a or newer with the **Optimization Toolbox** (Home > Add-Ons > Get Add-Ons); the
   official field optimization and integer rounding call `lsqnonlin`.
2. **Native build and Python environment** (do not use an Intel/Rosetta Homebrew; its libraries will
   not link). [environment.yml](environment.yml) lists the packages:
   ```bash
   conda env create -f environment.yml        # existing env: conda env update -n scrapshell -f environment.yml
   git submodule update --init --depth 1 third_party/Directional
   scripts/build_mesher.sh
   ```

## Run

```bash
scripts/run_all.sh                 # all eight models -> web/data/*.json (about 30 s per model)
scripts/run_all.sh metro vase      # a subset
scripts/serve.sh                   # regenerate web/index.html and serve it at http://localhost:8765/
scripts/snap.sh /tmp/shots "metro:4 metro:8"   # headless-Chrome screenshots of model:stage
```

From MATLAB: `addpath matlab; run_with_retries('metro')` or `run_phcpf('metro', 'upto', 4)`.
Intermediate files (raw field, primal / dual / planar OFF meshes, logs) go to `results/<model>/`.

## Layout

- `third_party/PH-CPF/`: official MATLAB code (CC BY-SA 4.0), without its Windows binaries.
  Local edits are listed in `third_party/PH-CPF/PATCHES.md`.
- `third_party/Directional/`: submodule, pinned to `2573873`.
- `tools/phcpf_mesher/`: replaces the official MeshCutter / Parametrization / Mesher executables.
- `matlab/`: driver (`run_phcpf.m`, `run_with_retries.m`), rounding glue, dualizer, planarizer,
  metrics, JSON export. `run_phcpf.m` also takes a per-face η and writes variant runs beside a
  model's baseline (`'out'`); `run_eta_checks.m` runs the extension's sizing checks on Metro.
- `web/index.html`: the viewer (generated by `scripts/build_web.sh` from `web/viewer.html`, the
  source, which is also what gets published as a claude.ai artifact); data in `web/data/`.
- `scripts/`: build, batch run, local server, screenshots.
- `stock/`: the stock of offcuts for the extension (Python). One deposit of 4000 scraps in meters,
  shared by every model, in `data/stock/deposit_s0.json`; every model is scaled to a 15 m span to
  use it. Shape descriptors (`describe.py`), four shape families (`shapes.py`), size distributions
  (`sizes.py`), model surfaces and baseline panels from `web/data` (`model_ref.py`), the generator
  (`python -m stock.generate`), the placement on a model (`python -m stock.place`) and the figures
  plus a quick fit preview against the baseline panels (`python -m stock.plot`, written to
  `results/stock/`).
- `sizing/`: stock-driven sizing for the extension (Python, in progress). `footprint.py` gives the
  panel PH-CPF aims for on each face, for one η or one per face, recovered from the viewer data
  (`python -m sizing.footprint`); `eta.py` writes per-face η files for `run_phcpf.m`
  (`python -m sizing.eta`), optionally from a smoothed size field (`--smooth`).
- `matching/`: panels against scraps (Python). `fit.py` is the exact fit test (any translation, 36
  rotations) with cheap filters; `assign.py` assigns panels to scraps, one each, and reports the fit
  rate and utilization (`python -m matching.assign metro`, written to `results/matching/`).

Tests are in `tests/`; run them with `~/miniconda3/envs/scrapshell/bin/python -m pytest` (or with
the env's `python` on your PATH).

## Differences from the official implementation

The official release only ships Windows executables for everything after the field optimization.
Those executables are replaced; every MATLAB step of the official pipeline runs unmodified.

**Unchanged official code:** curvature, guiding field, sizing, the CPF energies and their
Levenberg-Marquardt solve (`lsqnonlin`), the matching and singularities
(`compute_matching_and_sing`, passed to the cutter exactly as the official pipeline did), and the
integer rounding `BarrierPoissonSolve.m`. Parameters are those of `runme.m`. Code-vs-paper
discrepancies in the official code (sizing constant, an injectivity Jacobian typo, Eq. 16/27) are kept
as-is; see `third_party/PH-CPF/PATCHES.md`.

**Replaced executables** (same inputs and outputs, different implementation):

- *MeshCutter + Parametrization* -> `phcpf_mesher setup`: cut, comb and the seamless Poisson system
  from the current Directional (the official binaries used an older Directional). The system is
  exported in the variables `BarrierPoissonSolve.m` loads from `poisson.mat`. In the current
  Directional singular vertices are separate variables, so they are rounded together with the seam
  translations (a +-1/6 cone must sit on a lattice point anyway). `lengthRatio = 1/bbd`: the current
  Directional normalizes the field and the parameter length by the same factor, which reproduces the
  official intent of keeping the CPF's own magnitudes as the cell size.
- *Mesher* (CGAL) -> `phcpf_mesher mesh` (Directional's exact isoline mesher).
- *Dualizer* -> `matlab/dualize_bary.m`: barycentric dual, one boundary layer stripped (`CLIPPED_OLD`).
- *Planarizer* (Ceres) -> `matlab/planarize_hex.m`: supplement Eq. 23 with the paper's weights and the
  Eq. 14 barrier (s = 0.1 x previous length), 250 Levenberg-Marquardt iterations per pass, at most
  100 passes, lambda_p doubled per pass (capped at 1e6 as a numerical guard). Lambda_l is doubled per
  pass rather than every 5 solver iterations, since LM cannot change weights mid-solve. The unit
  normals use only their direction in the planarity term (`<n/|n|, e>`), which is what the
  supplement's homogeneous parametrization achieves in Ceres; with `<n, e>` plus a fixed-weight
  `|n| = 1` penalty, the growing lambda_p made collapsing a hard panel's normal to zero cheaper than
  flattening it, so a few panels stalled at 5-7% (and 0/0 renormalization caused NaN restarts).

**Robustness fixes** (needed on this platform, no change to the method):

- Directional's mesher merges coincident arrangement edges only when they have the same orientation
  (`generate_mesh.h`, the second clause of the test repeats the first). An isoline running exactly
  along a mesh edge, which happens whenever two adjacent singularities land on lattice points sharing
  a coordinate, then yields an inconsistent DCEL and an empty mesh. Fixed at build time in a local copy
  of the header (`tools/phcpf_mesher/CMakeLists.txt`); the submodule is untouched.
- The mesher can loop forever on a degenerate parameterization; each call is capped at 120 s.
- A failed rounding is retried with the paper's continuation (continuity weight doubled, up to 10
  times, Sec. 4.5); the official `runme.m` retried at most 3 times.
- `OPT.m`: removed the `CheckGradients` option that MATLAB R2026a rejects; `fast_sparse` is replaced
  by `sparse` (the official one is a Windows MEX file).

## Results

Six of eight models reach the paper's 1% planarity target; car and three holes stop at integer
rounding (the official `BarrierPoissonSolve` fails at every continuation level). Regenerate the table
with `python scripts/results_table.py` after `scripts/run_all.sh`; the same numbers are in the viewer.

<!-- results-table -->

| Model | Faces in | Doublings | LM iters | Singularities | Panels (hex + other) | Planarity max before → after | Planarity mean | Distance max / bbd |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Metro | 1418 | 1 | 129 | 21 | 422 + 14 | 4.60% → 0.87% | 0.122% | 0.0072 |
| Intersection | 1423 | 0 | 27 | 22 | 502 + 20 | 8.88% → 0.95% | 0.120% | 0.0050 |
| Six | 1365 | 0 | 55 | 5 | 412 + 4 | 5.60% → 0.77% | 0.171% | 0.0055 |
| Vase | 2404 | 0 | 31 | 0 | 394 + 0 | 2.14% → 0.85% | 0.251% | 0.0017 |
| Three holes | 1508 | 10 | 382 | — | — | stopped at stage 4 | — | — |
| Car | 1496 | 10 | 580 | — | — | stopped at stage 4 | — | — |
| Sphere with cones | 1972 | 0 | 62 | 12 | 702 + 12 | 2.72% → 0.86% | 0.252% | 0.0029 |
| C | 996 | 0 | 28 | 14 | 877 + 12 | 5.75% → 0.85% | 0.180% | 0.0062 |

Generated 2026-11-11 20:37; field solved with MATLAB lsqnonlin (Levenberg-Marquardt), as in the official code. Paper (Fig. 24/25): 255–1449 faces, planarity max 0.3–1.0% (target 1%), mean 0.03–0.49%, distance 0.001–0.015 × bbd.

## Credits

PH-CPF method and MATLAB code: Kacper Pluta, Michal Edelstein, Amir Vaxman, Mirela Ben-Chen
(github.com/michaled/PH-CPF, CC BY-SA 4.0). Directional: Amir Vaxman et al. (MPL 2.0).
