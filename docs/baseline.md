# Baseline system: PH-CPF on macOS

The baseline turns a triangle mesh of a target surface into a mesh of planar hexagonal panels that
stays close to the surface. It reproduces PH-CPF, *Planar Hexagonal Meshing using Coordinate Power
Fields* (Pluta, Edelstein, Vaxman, Ben-Chen, ACM TOG 40(4), 2021). The official MATLAB code runs
unmodified up to integer rounding, and our own code replaces the Windows-only executables after it.

This document explains how the system works, stage by stage, and where each piece lives. Setup and
run commands are in the [README](../README.md#setup-macos-apple-silicon). Current results are in
[progress.md](progress.md).

## The idea

A hexagonal panel mesh is the dual of a triangle mesh in which most vertices have six neighbours.
PH-CPF designs that triangle mesh directly on the surface:

1. On every triangle of the input mesh, two tangent vectors **U** and **V** map a regular triangular
   lattice in the plane onto the surface. They set the direction and spacing of the grid there.
2. If neighbouring triangles agree on the lattice, up to a rotation by a multiple of 60°, the field
   can be integrated into a seamless parameterization. The integer isolines of its three coordinate
   functions then cut the surface into triangles.
3. The dual of that triangle mesh is a hexagonal mesh: each vertex becomes a panel.

The field is designed so the hexagons come out almost planar. U and V are kept *conjugate* with
respect to the surface's curvature, and *sized* so each panel stays within a distance δ of the
surface. A final optimization then makes every panel planar to within 1%.

Curvature decides the panel shape. Dome-like (elliptic) regions get convex hexagons, saddle-like
(hyperbolic) regions get non-convex bow-tie hexagons, and cylinder-like (parabolic) regions get
elongated brick hexagons.

## Pipeline at a glance

The stage numbers match the viewer and `run_phcpf(id, 'upto', k)`.

| # | Stage | Output | Code | Origin |
| ---: | --- | --- | --- | --- |
| 1 | Input | Triangle mesh | `third_party/PH-CPF/data/*.off` | Official data |
| 2 | Curvature regions | Elliptic, parabolic, hyperbolic or planar, per face | `SHAPEOP.guiding_field2` | Official MATLAB |
| 3 | Guiding field and sizing | Guiding direction field; sizing metric `g` | `CPF.set_param`, `SHAPEOP.curvature_sizing` | Official MATLAB |
| 4 | Coordinate power field | U, V per face; matching; singularities | `CPF.opt_power_cfs`, `CPF.compute_matching_and_sing` | Official MATLAB |
| 5 | Seamless parameterization | Coordinate functions h1..h3, rounded to integers across cuts | `phcpf_mesher setup` + `BarrierPoissonSolve.m` | Ours (C++) + official MATLAB |
| 6 | Primal mesh | Triangle mesh from integer isolines | `phcpf_mesher mesh` | Ours (C++, Directional) |
| 7 | Hexagonal dual | Polygon mesh, mostly hexagons | [matlab/dualize_bary.m](../matlab/dualize_bary.m) | Ours (MATLAB) |
| 8 | Planarization | Planar panels near the surface | [matlab/planarize_hex.m](../matlab/planarize_hex.m) | Ours (MATLAB) |

[matlab/run_phcpf.m](../matlab/run_phcpf.m) drives all eight stages for one model.
[matlab/run_with_retries.m](../matlab/run_with_retries.m) wraps it with the continuation described in
stage 5.

## Parameters

The defaults are those of the official `runme.m`.

| Parameter | Value | Meaning |
| --- | --- | --- |
| `n` | 6 | Symmetry of the field: 6 gives hexagons, 4 gives quads |
| `delta` (δ) | 0.005 | Maximum distance of a panel from the surface, as a fraction of the bounding-box diagonal (bbd) |
| `eta` (η) | 0.002 | Maximum panel area, as a fraction of the total surface area |
| `mp` | 1 | Planarity target, in percent (paper Eq. 29) |
| `lsm`, `lc`, `lsi`, `la`, `lin`, `lcj`, `lo` | 1 each | Weights of the field energies (stage 4) |
| `siz`, `optinit` | `'k'` | Curvature sizing; field initialized from the guiding field |
| `iter_opt` | 0..10 | Continuity-weight doublings (stage 5 continuation) |

## Stage by stage

### 1. Input

An OFF triangle mesh from `third_party/PH-CPF/data/`. The eight models are listed in
[matlab/model_catalog.m](../matlab/model_catalog.m): four open shells (metro, intersection, six_s0,
vase) and four closed models (3holes_r, car_r, sphere2con_s0, c_s0). The bounding-box diagonal (bbd)
sets the length scale for δ and for the distance metric.

### 2. Curvature regions

The official shape operator gives, per face, the principal curvatures κ_min and κ_max and their
directions. `SHAPEOP.guiding_field2` summarizes them as a magnitude ρ and an angle φ that encodes
the curvature ratio. φ = π/4 is parabolic, smaller is hyperbolic, larger is elliptic. We classify
faces with the official thresholds: planar if ρ ≤ 0.01, parabolic if |φ − π/4| < 0.05, then
hyperbolic or elliptic. This stage is for display only. The field uses the curvature directly.

### 3. Guiding field and sizing metric

**Guiding field.** A smooth direction field, mostly aligned with the direction of least absolute
curvature. It initializes the frame field and the alignment energy pulls toward it.

**Sizing metric.** `SHAPEOP.curvature_sizing` builds a 2×2 metric $g_f$ on every face $f$
([SHAPEOP.m:259](../third_party/PH-CPF/code/SHAPEOP.m#L259)). In the principal frame,

$$g_f = \mathrm{diag}\left(\tilde\kappa_\min / \delta,\ \tilde\kappa_\max / \delta\right), \qquad \tilde\kappa = \sqrt{\kappa^2 + \varepsilon^2}, \quad \varepsilon = \delta / \eta .$$

The ellipse $\{u : u^\top g_f u = 1\}$ has semi-axes $\sqrt{\delta / \tilde\kappa}$ along the principal
directions. At that distance from the tangent point, the surface has moved away from its tangent plane by
$\tfrac12 \kappa x^2 = \delta/2$. So the ellipse is the region a flat panel can cover while staying
within δ/2 of the surface: the panel's **footprint**. The smoothing ε caps the footprint where a
principal curvature vanishes, in flat directions and along the asymptotic directions of saddles,
where the panel could otherwise grow without bound. There, the size falls back to $\sqrt{\eta}$, the
maximum-area cap.

The viewer draws these ellipses in stage 3. This metric is the hook for the stock-driven extension;
see [plan.md](plan.md#key-observation-the-footprint-is-the-sizing-metric).

### 4. Coordinate power field

**Unknowns.** Per face, the 2D tangent vectors U and V in the face's local basis. Per interior edge,
a complex number z: the edge expressed in lattice coordinates.

**Energies.** Levenberg–Marquardt (`lsqnonlin`) minimizes a weighted sum of least-squares terms
([CPF.m:649](../third_party/PH-CPF/code/CPF.m#L649) and after):

| Term | Weight | What it asks for |
| --- | --- | --- |
| Smoothness | `lsm` | U and V vary smoothly between neighbouring faces (supplement Eq. 19) |
| Continuity | `lc` | On each interior edge, both faces agree on the edge's lattice coordinates z, compared through their n-th powers. Raising to the n-th power removes the 60° ambiguity, hence "power field". This is the integrability condition |
| Sizing | `lsi` | $U^\top g U = V^\top g V = 1$: the lattice fills the footprint ellipse ([CPF.m:901](../third_party/PH-CPF/code/CPF.m#L901)) |
| Conjugacy | `lcj` | $U^\top S V = 0$ for the shape operator S: the condition that makes the dual hexagons planarizable |
| Orthogonality | `lo` | U · V = 0 |
| Alignment | `la` | U follows the guiding field and the boundary |
| Injectivity | `lin` | Barrier that keeps det[U V] > 0, so no lattice cell flips |

The field typically converges in 27–129 iterations (up to 580 on the failing models).

**Matching and singularities.** `compute_matching_and_sing` finds, on every interior edge, the
multiple of 60° that best aligns the two faces' lattices: the *matching*. It then finds the vertices
around which the matchings do not add up to zero: the *singularities*. Singular vertices become
pentagons or heptagons in the final mesh. These are passed unchanged to the cutter, as in the
official pipeline.

### 5. Seamless parameterization and integer rounding

1. **Export from MATLAB.** The field becomes the gradients of six "unit-root" functions per face
   (`CPF.grads2rawfield`), written to `<id>.rawfield`. The matching and singularities go to
   `<id>_matching.txt` and `<id>_singin.txt`.
2. **Cut and set up** (`phcpf_mesher setup`, [tools/phcpf_mesher/main.cpp](../tools/phcpf_mesher/main.cpp)).
   The mesh is cut open along a graph through the singularities, the field is combed, and the
   seamless Poisson problem is assembled: find coordinate functions whose gradients match the field,
   with offsets across each cut that must become integers. The system is exported in the variables
   the official `BarrierPoissonSolve.m` expects (A, b, C, G, x0, integer indices, the cut mesh).
   Parameter length ratio = 1/bbd, which keeps the CPF's own magnitudes as the cell size.
3. **Round** ([matlab/barrier_rounding.m](../matlab/barrier_rounding.m) runs the official
   `BarrierPoissonSolve.m` unmodified). Integer variables are rounded one at a time. After each, the
   rest is re-solved with Levenberg–Marquardt under a barrier that keeps every face's six gradient
   vectors in their cyclic order, so no triangle flips. The integer variables are the seam offsets
   and, in the current Directional, the singular vertices' values: a ±1/6 cone must sit on a lattice
   point.
4. **Continuation.** If rounding fails, `run_with_retries` reruns the field with `iter_opt + 1`,
   which doubles the continuity weight during warm-up passes and makes the field closer to integrable
   (paper Sec. 4.5). It retries up to 10 times. Car and Three holes fail at every level.

### 6. Primal mesh from isolines

`phcpf_mesher mesh` takes the rounded solution and runs Directional's isoline mesher. In every
triangle it intersects the integer isolines of the three independent coordinate functions, using
exact rational arithmetic (GMP). It then stitches the pieces across the cuts into one mesh. On the
closed models every face is a triangle. On the open shells, 7–9% of faces are quads (and the odd
pentagon): triangles clipped by the surface boundary.

Two safeguards:

- A **patch** to Directional's `generate_mesh.h`, applied to a build-local copy in
  [tools/phcpf_mesher/CMakeLists.txt](../tools/phcpf_mesher/CMakeLists.txt). Coincident edges were
  merged only when they had the same orientation (the test repeated one clause). An isoline running
  exactly along a mesh edge in the opposite direction then produced an inconsistent half-edge
  structure and an empty mesh. This happens whenever two adjacent singularities round to lattice
  points sharing a coordinate.
- A **120 s cap** per call, because the mesher can loop forever on a degenerate parameterization.

### 7. Hexagonal dual

[matlab/dualize_bary.m](../matlab/dualize_bary.m) makes every interior primal vertex a panel through
the centroids of its incident triangles, ordered around the vertex and oriented like the surface.
Valence-6 vertices give hexagons; singular vertices give pentagons or heptagons. Panels whose
one-ring touches the boundary are dropped (`strip = 1`, like the official dualizer's clipped mode),
because they are built from the clipped boundary faces and come out irregular.

Before planarization, the panels are 2–9% non-planar (Eq. 29 maximum).

### 8. Planarization

[matlab/planarize_hex.m](../matlab/planarize_hex.m) follows supplement Eq. 23. The unknowns are the
vertex positions and one normal per panel. Coordinates are rescaled to unit mean edge length during
the solve.

| Term | Weight | What it asks for |
| --- | --- | --- |
| Planarity | λp = 0.01, doubled every pass (capped at 1e6) | Every panel edge orthogonal to the panel normal. Only the normal's direction counts |
| Unit normal | λn = 1 | Gauge term for the normal's length |
| Distance | λd = 0.5 | Each vertex stays near its closest point on the input surface (recomputed every pass) |
| Central symmetry | λs = 1 | Opposite vertices of even-sided panels stay symmetric about the panel centre |
| Normal drift | λc = 0.01 | Vertices move little along their vertex normal (averaged from the panel normals) relative to the previous pass |
| Length barrier | λl = 0.001, doubled every pass | Edges and diagonals do not collapse (paper Eq. 14, s = 0.1 × previous length) |

Each pass runs up to 250 Levenberg–Marquardt iterations. Passes repeat, with the weights doubled,
until the worst panel is planar to within `mp` = 1%, or for at most 100 passes. The successful
models take 11–16 passes.

## Metrics

Computed by [matlab/metrics_phcpf.m](../matlab/metrics_phcpf.m) and
[matlab/planarity_eq29.m](../matlab/planarity_eq29.m).

- **Planarity** (paper Eq. 29), in percent. For every run of four consecutive vertices v1..v4 of a
  panel, take n as the unit normal of the two diagonals, (v3 − v1) × (v4 − v2). The error is
  |⟨n, v2 − v1⟩| divided by the mean diagonal length. A panel's value is 100 × the RMS over its runs.
  We report the maximum (the target is ≤ 1%) and the mean.
- **Distance.** One-sided distance from panel vertices and panel centroids to the input surface,
  divided by bbd. We report the maximum and the mean.
- **Counts.** Panels, hexagons versus other polygons, singularities, iterations, time per stage.

## Files and data flow

### Code

| Path | Role |
| --- | --- |
| [matlab/run_phcpf.m](../matlab/run_phcpf.m) | Driver: runs stages 1–8 and exports everything for the viewer |
| [matlab/run_with_retries.m](../matlab/run_with_retries.m) | Continuation: reruns with `iter_opt` = 1..10 until stage 8 is reached |
| [matlab/barrier_rounding.m](../matlab/barrier_rounding.m) | Loads the exported Poisson system and runs the official `BarrierPoissonSolve.m` |
| [matlab/dualize_bary.m](../matlab/dualize_bary.m) | Stage 7 |
| [matlab/planarize_hex.m](../matlab/planarize_hex.m) | Stage 8 |
| [matlab/closest_point_mesh.m](../matlab/closest_point_mesh.m) | Closest points on the input mesh (planarizer and distance metric) |
| [matlab/metrics_phcpf.m](../matlab/metrics_phcpf.m), [matlab/planarity_eq29.m](../matlab/planarity_eq29.m) | Metrics |
| [matlab/phcpf_recorder.m](../matlab/phcpf_recorder.m) | Records every Levenberg–Marquardt iterate of the field for the viewer |
| [tools/phcpf_mesher/](../tools/phcpf_mesher/) | C++ replacement for MeshCutter, Parametrization and Mesher |
| [third_party/PH-CPF/](../third_party/PH-CPF/) | Official MATLAB code. Local edits are in [PATCHES.md](../third_party/PH-CPF/PATCHES.md) |
| `third_party/Directional/` | Submodule, pinned to `2573873` |
| [scripts/](../scripts/) | Build, batch run, local server, screenshots, results table |

### Intermediate files

Each run writes to `results/<id>/` (git-ignored):

| File | Contents |
| --- | --- |
| `<id>.rawfield` | Field as six gradient vectors per face (stage 4 → 5) |
| `<id>_matching.txt`, `<id>_singin.txt` | MATLAB's matching and singularities, passed to the cutter |
| `<id>_A.txt`, `_b.txt`, `_C.txt`, `_G.txt`, `_x0.txt`, `_int.txt`, `_singidx.txt`, `_x2c.txt`, `_cutV.txt`, `_cutF.txt`, `_FN.txt`, `_meta.json` | Seamless Poisson system exported by `phcpf_mesher setup` |
| `<id>_poisson.mat` | The same system in the format `BarrierPoissonSolve.m` loads |
| `<id>_x.txt` | Rounded solution |
| `<id>_corner.txt`, `<id>_sing.txt`, `<id>_stats.json` | Coordinate functions per triangle corner, singularities used, mesher statistics |
| `<id>_tri.off`, `<id>_dual.off`, `<id>_planar.off` | Primal mesh, hexagonal dual, planarized panels |

Logs go to `results/<id>.log`.

### Viewer data

`web/data/<id>.json` holds the geometry of every stage, including up to 36 animation frames each for
the field optimization and the planarization passes. `web/data/<id>.meta.json` holds the metrics and,
for a failed run, the error and the last stage reached. `web/data/manifest.json` lists the models.

## Viewer

[web/viewer.html](../web/viewer.html) is the source. [scripts/build_web.sh](../scripts/build_web.sh)
wraps it into `web/index.html`, which GitHub Pages serves at
[zincr.github.io/ScrapShell/web](https://zincr.github.io/ScrapShell/web/). The empty `.nojekyll` file
stops Pages from running Jekyll over the repository. Run it locally with `scripts/serve.sh`. Take
headless screenshots with `scripts/snap.sh OUTDIR "metro:4 metro:8"`.

## How it differs from the official release

The [README](../README.md#differences-from-the-official-implementation) and
[PATCHES.md](../third_party/PH-CPF/PATCHES.md) have the full list. In short:

- **Unchanged official code:** curvature, guiding field, sizing, the field energies and their
  solve, matching and singularities, and `BarrierPoissonSolve.m`, all with the `runme.m` parameters.
  Known discrepancies between the code and the paper (sizing constant, an injectivity typo, Eq. 16
  and 27) are kept as-is for fidelity.
- **Replaced executables:** MeshCutter and Parametrization → `phcpf_mesher setup` (current
  Directional). Mesher → `phcpf_mesher mesh`. Dualizer → `dualize_bary.m`. Planarizer (Ceres) →
  `planarize_hex.m`. The planarizer doubles λl per pass instead of every 5 solver iterations,
  because Levenberg–Marquardt cannot change weights mid-solve.
- **Fixes:** the mesher's reversed-edge merge; normal direction only in the planarity term; a 120 s
  mesher cap; up to 10 continuation levels instead of 3; `OPT.m` drops an option MATLAB R2026a
  rejects; `fast_sparse` (Windows MEX) → `sparse`.

## Limitations that matter for the extension

- **One size rule everywhere.** Sizing is an equality, $U^\top g U = 1$, so every panel is as large
  as its footprint allows, capped by η where the surface is flat. The only controls are the global
  δ and η, so nothing about the available stock can influence panel size.
- **Panel shape is fixed by curvature.** There is no notion of material or stock anywhere.
- **Closed models with many singularities fail at rounding** (Car, Three holes). Anything that makes
  the field less integrable, such as rapidly varying panel sizes, risks the same failure.
- **Planarization dominates runtime** on the larger models (55–66% on Metro, Intersection, Six and C).
  A full run still takes 7–20 s, so iterating the pipeline in a loop is affordable.
