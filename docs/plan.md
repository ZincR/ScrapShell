# Project plan: Scrap Shells

Proposed 2026-11-21. Revised 2026-11-24: a per-face cap η replaces the size factor σ (see
[Why the cap varies per face](#why-the-cap-varies-per-face)). The final due date is not yet
confirmed (see [Open decisions](#open-decisions)).
For where things stand, see [progress.md](progress.md). For how the baseline works, see
[baseline.md](baseline.md).

## Question

Can a stock of irregular flat offcuts be assembled into a shell that approximates a target surface?
We start from a planar-panel meshing method (PH-CPF) and adapt it so that the stock informs the panel
layout: big scraps make big panels where the surface allows, small scraps make small ones.

## Objectives

From the proposal. A layout succeeds when it meets the three core objectives.

1. **Fit:** every panel fits inside the offcut it is assigned to.
2. **Geometry:** panels approximate the target surface in distance and follow its curvature.
3. **No gaps:** neighbouring panels share edges and form a closed, manifold surface.

If time allows:

4. **Fewest pieces:** use as few offcuts as possible.
5. **Least waste:** leave as little of each used offcut as possible.

## Approach

### Key observation: the footprint is the sizing metric

The proposal's working hypothesis is "the surface sets the footprint, the stock sets the fill". A
footprint is the region around a point where the tangent plane stays within ε of the surface.

PH-CPF already computes exactly this. Its per-face sizing metric g describes an ellipse with
semi-axes $\sqrt{\delta/\tilde\kappa}$ along the principal directions. Over that ellipse the tangent
plane stays within δ/2 of the surface (see [baseline.md, stage 3](baseline.md#3-guiding-field-and-sizing-metric)).
PH-CPF even caps the footprint where a curvature vanishes, which the proposal also called for. The
field energy then asks every panel to fill its footprint exactly: $U^\top g U = 1$.

So the extension needs only one new input to the field: the cap, set per face from the stock.
PH-CPF uses one cap for the whole surface, $\epsilon = \delta/\eta$, where η is the maximum element
area. We give each face its own $\eta_f$, so its footprint has semi-axes
$a_i = \sqrt{\delta/\sqrt{\kappa_i^2 + \epsilon_f^2}}$ with $\epsilon_f = \delta/\eta_f$. δ, and with
it the distance to the surface, does not change. Where the surface is flat, $\eta_f$ sets the
panel's size; where it is parabolic, $\eta_f$ clips the panel's long axis. A panel can be smaller than
its curvature allows, never larger.

### Why the cap varies per face

The first draft of this plan (2026-11-21) kept PH-CPF's single η and scaled g by a per-face factor
σ ≤ 1 instead. Measuring the footprints of the baseline runs (`python -m sizing.footprint`,
[sizing/footprint.py](../sizing/footprint.py)) ruled that out on 2026-11-24:

- **The cap, not the curvature, sets the panel size almost everywhere.** With the baseline's η it
  decides 65–100% of each finished model's surface (Metro 98%). σ ≤ 1 can only shrink panels, so
  panels could not grow even where the scraps are large.
- **Raising η for every face stretches the panels.** That is the only way to let σ-scaled panels
  grow, but a larger η uncaps the flat direction first. The 90th-percentile aspect (long ÷ short
  semi-axis) rises from 1.1–2.2 at the baseline's η to 2.0–5.7 at 16× (Metro 1.7 to 4.3). σ scales
  both axes alike, so it cannot undo the stretch.
- **Stretched panels waste the stock.** Taking the largest panel each scrap holds at a given aspect,
  sliding and turning it (the capacity test of item 2), the deposit supplies panels covering 50% of its area at aspects 1.5–2,
  42% for round panels and 31% at 5:1. At 5:1 only 2% of scraps hold a Metro-sized panel (0.275 m²),
  against 10% at 1.5–2. Panels and scraps match best when their shapes match: 3.5–5:1 boards use
  53% of their area with a 5:1 panel and 16% with a round one.

A per-face η lets each face grow where the stock allows, and sets its stretch anywhere between
round and what its curvature allows. One number sets both, so the size field also keeps the aspect
near 1.5–2 (item 2).

### Pipeline

| Proposal step | Planned implementation | Status |
| --- | --- | --- |
| **Seed:** one footprint per scrap, sized by the stock | Size field η over the surface (the per-face cap), computed by 1D optimal transport between the deposit's scrap capacities and the surface's room for panels | Footprint model done; rest to do |
| **Grow:** regions grow toward their scrap's size | PH-CPF field and seamless parameterization with the per-face cap | Baseline done; η hook to do |
| **Intersect:** flat panels that close the shell | PH-CPF dual and planarization. Closed by construction | Done (baseline) |
| **Match:** pair each scrap with a panel it fits | Exact fit test plus linear assignment between panels and scraps | To do |
| **Check / Revise:** reseed where it fails | Feedback loop: lower η where panels fit no scrap | To do |
| **Cut:** trim each panel from its scrap | Place each panel in its scrap; export the cut outlines | To do |

### Why this instead of the proposal's pipeline

The proposal planned power diagrams for seeding, centroidal Voronoi cells grown on geodesics, and
tangent plane intersection (Zimmer et al. 2012) to close the shell. We change course for four reasons.

1. **It builds on the baseline.** The proposal's pipeline would be a second, separate system that
   uses none of the PH-CPF reproduction.
2. **"No gaps" comes for free.** The hexagonal dual of the isoline triangle mesh is a closed polygon
   mesh. The proposal named coverage as its main risk; this route removes it.
3. **Optimal transport keeps its role.** The size field is a 1D transport problem, and Match is a
   linear assignment, a discrete transport problem.
4. **It fits the remaining time.** A full pipeline run takes 7–20 s, so a fit-and-resize loop costs
   minutes.

Tangent plane intersection stays in the related work, and as a stretch comparison.

## Work plan

### Phase 0: baseline (done)

PH-CPF on macOS, 6 of 8 models at 1% planarity, web viewer. See [progress.md](progress.md).

### Phase 1: stock and stock-driven sizing

1. **Stock generator** (done, [stock/](../stock/); details in [stock.md](stock.md)). One deposit
   shared by every model, in meters.
   - [data/stock/deposit_s0.json](../data/stock/deposit_s0.json): 4000 scraps, exactly 1000 from each
     family (remnant, shard, blob, board), each family with the same size distribution. Size
     (√area) is lognormal with median 0.5 m, cut off at 0.15–1.5 m: 1190 m² in all.
   - Every model is given the same span: its bounding-box diagonal is scaled to 15 m. δ stays 0.5%
     of the diagonal, which is then 75 mm on every model. At that span the models' surfaces are
     75–284 m², so the deposit holds 4–16× each model's area and 4.5–10 scraps per baseline panel:
     piece count never limits a layout.
   - **Fit, not count, is the constraint, but only for panels larger than the baseline's.** A
     quick fit preview found that only 7–38% of scraps can hold a typical baseline panel centred on
     their incenter (Metro 7%). The exact test (item 5) finds far more: every one of Metro's 436
     baseline panels fits a scrap of its own, and fit starts to bind between a 0.5 and a 1 m² median
     panel (item 6). We keep the 0.5 m median: the size field still has work to do, choosing how
     large panels can grow.
   - `python -m stock.generate` writes the deposit; `python -m stock.place` prints each model's
     physical size against it; `python -m stock.plot` draws the figures and runs the fit preview.
2. **Size field η by 1D optimal transport.** Python, in [sizing/](../sizing/).
   - **Footprint model** (done, [sizing/footprint.py](../sizing/footprint.py)). It recovers each
     face's principal curvatures, in 1/m, from the sizing ellipses in the viewer data, so no MATLAB
     run is needed, and gives the footprint for any η: one value, or one per face. A panel is a
     regular hexagon in g with circumradius 0.70, fitted on the six finished baseline runs (0.68–0.71
     on each). This predicts each baseline's panel count within 6%.
   - **Scrap capacity.** For each scrap, the area of the largest panel it holds at a given aspect,
     anywhere inside it and at any rotation: squash the scrap by the aspect along the panel's long
     axis and take its largest inscribed circle. Sliding is exact; rotation is sampled. This replaces
     the first draft's inscribed-circle rule, which passed 1.6–3.4× too few scraps
     ([stock.md](stock.md#fit-preview)).
   - **Monotone coupling.** Sort faces by their room, the largest panel their curvature allows with an
     aspect of at most 2 (unbounded where the surface is flat), and scraps by capacity, both largest
     first. Walk the two lists together: each scrap claims the surface area of one panel of its size,
     and those faces take that panel, or their room if smaller, as their target. Big scraps go where
     the surface has room. The deposit holds far more than one surface needs, so this is a partial
     transport that also chooses which scraps to use. Check Hall's condition: for every size, at least
     as many scraps hold a panel of that size as the surface needs panels at least that large.
   - **From target to η** (done 2026-11-24, `eta_for_area` in [sizing/eta.py](../sizing/eta.py)).
     On each face, solve for $\epsilon_f$ so the panel has the target area (a quadratic in
     $\epsilon_f^2$), then raise $\epsilon_f$ where needed to keep the aspect at most 2 (also closed
     form). A target larger than the curvature allows stops at η = 16× the baseline's (3.2%), the
     top of the ramp that rounded.
   - **Smooth the target, not just η** (done 2026-11-24, `smooth_eta`). Smooth the log target size
     with Laplacian steps over the face adjacency, then convert. Because the conversion matches each
     face's panel area to its target, a smooth target gives a smooth metric even where the curvature
     estimate jumps from face to face. Metro's faces are small next to its panels (about 0.3 m
     against a 0.68 m panel radius), so it takes many steps: 400 bring the median log-size jump
     across an edge from 0.027 to 0.007, below the baseline's 0.008, at about the same panel count.
     A hard limit on the jump per edge, as first planned, backfires: the few faces where the
     curvature allows only small panels pull the whole field down in cones around them (620–950
     panels instead of 280), so it is left out.
   - *Done when* η can be plotted on metro.
3. **Hook η into CPF** (done 2026-11-24).
   - `CPF.set_param` takes η as one value or one per face, each a fraction of the total area. It
     checks the length and uses `cpf.delta ./ cpf.eta` for `k_eps`; `SHAPEOP.curvature_sizing` and
     `check_curvature_sizing` use `e.^2` for `e^2`
     ([PATCHES.md](../third_party/PH-CPF/PATCHES.md)). g is built inside `set_param`, so `mel2d`,
     the weights it normalizes and `opt_init` all follow.
   - [run_phcpf.m](../matlab/run_phcpf.m) takes `'eta'` as one value, one per face, or a text file
     with one per face ([sizing/eta.py](../sizing/eta.py) writes them), and stores a per-face η in
     the viewer data as a packed array. Its `'out'`, `'label'` and `'group'` options write a variant
     run beside the model's baseline; the manifest lists variants after the catalog's models, and
     `stock.place.model_ids` leaves them out unless asked.
4. **Sanity checks on metro** (run 2026-11-24, [run_eta_checks.m](../matlab/run_eta_checks.m)).
   Every run is in the viewer's "Sizing checks" group.

   | Run | η | Result | Continuation level | Singularities | Panels | Planarity max | Distance max ÷ bbd |
   | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
   | Baseline | 0.2% | Stage 8 | 1 | 21 | 436 | 0.87% | 0.0072 |
   | η per face = baseline | 0.2% per face | Stage 8, identical to the baseline | 1 | 21 | 436 | 0.87% | 0.0072 |
   | η 4× | 0.8% per face | Rounding fails at every level | 0–10 | 41–54 | – | – | – |
   | η 1× to 16× | 0.2% to 3.1%, log-linear along the long axis | Stage 8 (one mesher crash at level 4) | 5 | 47 | 233 | 0.64% | 0.0076 |
   | η 4×, aspect ≤ 2 | ≤ 0.8% per face ([sizing/eta.py](../sizing/eta.py)) | Rounding fails at every level | 0–10 | 46–62 | – | – | – |
   | δ and η 4× | 0.8%, with δ 2% | Stage 8 at the first try | 0 | 6 | 89 | 0.65% | 0.0193 |
   | η 4×, smoothed | 0.29% to 3.2% per face (median 0.73%): 4×'s panel sizes after 400 smoothing steps ([sizing/eta.py](../sizing/eta.py)) | Stage 8 | 4 | 31 | 224 | 0.80% | 0.0095 |
   | η 4×, smoothed, aspect ≤ 2 | 0.065% to 3.2% per face (median 0.47%), the same with the aspect cap | Stage 8 | 8 | 29 | 261 | 0.49% | 0.0106 |

   - **The hook is exact.** Per-face η equal to the baseline's reproduces the baseline run bit for
     bit, and on every run MATLAB's sizing ellipses match [sizing/footprint.py](../sizing/footprint.py)
     to within 5 × 10⁻⁶.
   - **A larger η breaks rounding on Metro, varying or not.** A constant 4× never rounds; the ramp
     rounds only after five doublings of the continuity weight. Both carry more than twice the
     baseline's singularities.
   - **Not panel size:** growing δ with η (the same footprint shapes, twice as large) rounds at once,
     with fewer singularities than the baseline. **Not stretch alone:** capping the aspect at 2 still
     fails.
   - **The cause is a rough metric, and smoothing the size field fixes it.** At the baseline's η the
     cap evens out the footprints. With a larger η their size follows the per-face curvature
     estimate, which jumps between neighbouring faces: the median change in log footprint size
     across an edge is 0.027 at 4× and 0.026 with the aspect cap, against 0.008 at the baseline and
     with δ and η grown together. Smoothing 4×'s panel sizes (item 2) brings it to 0.007 at about the
     same panel count, and that run rounds after four doublings, with 31 singularities instead of
     41–54. With the aspect cap the jump is 0.010 (where the cap binds, the size follows the
     curvature again) and rounding needs eight doublings.
   - **Smoothed, larger panels keep the geometry.** Both smoothed runs reach the 1% planarity target
     (0.80% and 0.49%) and stay within 0.011 bbd of the surface, between the baseline (0.007) and
     δ grown with η (0.019).
   - **Still costly:** a size field that departs from the baseline's needs several doublings of the
     continuity weight. The size field should therefore be smoothed before every run; growing δ with
     η remains the fallback where rounding still fails.

### Phase 2: matching, cutting and feedback

5. **Fit test** (done 2026-11-24, [matching/fit.py](../matching/fit.py), with shapely).
   - Flatten each planarized panel P onto its best-fit plane. At ≤ 1% planarity this is safe.
   - P fits scrap S if, for some rotation, there is a translation that puts P inside S. That
     translation set is S shrunk by P: S minus (S's boundary ⊕ −P), a Minkowski sum built from the
     convex hulls of each boundary edge with −P. 36 rotations.
   - Filters first: a pair cannot fit if the scrap is smaller than the panel in area, inradius,
     circumradius or minimum width, and surely fits if the scrap's inscribed circle holds the
     panel's enclosing circle. On Metro's baseline this leaves 482k of 1.7M pairs for the exact
     test, which takes about 5 minutes on 12 cores.
   - Works for non-convex bow-tie panels too (split into triangles; 20% of Metro's panels).
   - Optionally grow every panel by the saw kerf before testing (not done yet).
6. **Assignment** (done 2026-11-24, [matching/assign.py](../matching/assign.py),
   `python -m matching.assign <run>`). `scipy.optimize.linear_sum_assignment`. The cost of a pair is
   the waste (scrap area minus panel area), with a cost for pairs that do not fit larger than any
   total waste, so it first fits as many panels as possible (checked against a maximum bipartite
   matching). Outputs, in `results/matching/<run>.json`: the assignment, the panels that fit
   nowhere, utilization, and each panel's placement in its scrap.

   First results, Metro at the deposit's 15 m span, one deposit of 4000 scraps:

   | Run | Panels | Median panel area | Fit rate | Utilization | Scraps that fit some panel | Inscribed-circle rule alone |
   | --- | ---: | ---: | ---: | ---: | ---: | ---: |
   | Baseline (η 0.2%) | 436 | 0.26 m² | 100% | 77% | 2579 | 45% |
   | η 1× to 16× | 233 | 0.42 m² | 100% | 70% | 3914 | 26% |
   | η 4×, smoothed, aspect ≤ 2 | 261 | 0.43 m² | 100% | 72% | 3146 | 15% |
   | η 4×, smoothed | 224 | 0.49 m² | 100% | 68% | 2383 | 17% |
   | δ and η 4× | 89 | 1.04 m² | 20% | 71% | 44 | 0% |

   - **At PH-CPF's own panel size every panel fits.** The first fit preview's 7% (a typical panel
     centred on the scrap's incenter) was far too pessimistic: with sliding and rotation 64% of the
     scraps hold some Metro panel. The inscribed-circle rule alone would fit only 45% of the panels.
   - **Panels can about double in area before fit binds.** Up to a 0.5 m² median every panel still
     fits, and utilization drops from 77% to 68–72%. At 1 m² only 44 scraps hold any panel.
   - **For the experiments (Phase 3):** Method A at the baseline η already fits every panel on
     Metro, so the comparison is how few pieces each method reaches at a 100% fit rate, and at what
     utilization, rather than the fit rate at one size.
7. **Feedback loop.** Lower η around panels that fit no scrap (or drop already-used scraps from the
   transport step) and rerun stages 3–8. Stop when every panel fits.
8. **Viewer stage 9.** Panels coloured by fit slack or scrap. A second view shows the pile of scraps
   with each panel's cut outline drawn inside its scrap.

### Phase 3: experiments

| Item | Choices |
| --- | --- |
| Models | metro, intersection (architectural shells), vase, sphere2con; six and c if time allows |
| Stock | The one deposit (4000 scraps in meters), with every model at a 15 m span. Optional ablation: a deposit with a different size spread |
| Method A | Plain PH-CPF with one global size. Sweep η and keep the best fit rate, then run the same assignment |
| Method B | Stock-driven sizing, one pass, then assignment |
| Method C | Stock-driven sizing with the feedback loop |

| Objective | Metric |
| --- | --- |
| Fit | % of panels assigned to a scrap they fit inside |
| Geometry | Planarity max and mean (Eq. 29); max and mean distance ÷ bbd; number of singularities |
| No gaps | Every interior edge shared by exactly two panels; boundary edges only on the surface boundary |
| Fewest pieces | Number of scraps used |
| Least waste | Utilization = total panel area ÷ total area of the scraps used |
| Cost | Runtime; number of feedback iterations |

**Expected headline:** A against C. One global size either wastes the big scraps or leaves panels
that fit nowhere. Stock-driven sizing should fit more panels with less waste, at the same planarity
and distance.

### Phase 4: write-up

- Final report: method, the footprint-as-metric observation, experiments, limitations.
- Presentation and figures: viewer screenshots via `scripts/snap.sh`, the pile view, the A-against-C
  comparison.
- Update [progress.md](progress.md) and the [README](../README.md).

## Timeline

Assumes about three weeks; adjust once the due date is confirmed.

| Week | Dates | Work |
| --- | --- | --- |
| 1 | Nov 21 – Nov 27 (Thanksgiving Nov 26) | Phase 1: stock generator, size field, η hook, sanity checks |
| 2 | Nov 28 – Dec 4 | Phase 2: fit test, assignment, feedback loop, viewer stage 9 |
| 3 | Dec 5 – Dec 11 | Phase 3 experiments; Phase 4 report and presentation |

The work splits cleanly between two people: Python (stock, size field, fit test, assignment: items
1, 2, 5 and 6) and MATLAB (η hook, sanity checks, feedback driver: items 3, 4 and 7).

## Risks

| Risk | Mitigation |
| --- | --- |
| Larger panels make the field harder to integrate: more singularities, rounding failures (the same failure as Car and Three holes). **Measured on Metro (item 4):** 4× η fails at every continuation level, with or without the aspect cap; the 1×–16× ramp passes after five doublings; 4×'s sizes smoothed pass after four (eight with the aspect cap) | Smooth the target size so the metric is smooth (item 2; done, and it makes 4× round); reuse the continuation in [run_with_retries.m](../matlab/run_with_retries.m); fall back to growing δ with η, which rounds at once; run experiments on the models that already succeed |
| The grid only roughly hits the target size, so first-pass panels miss their scraps | The feedback loop (item 7) |
| η sets a face's size and stretch together, so a parabolic face cannot have a large, round panel | Cap the aspect at 2, where the deposit fits panels best; if that limits too much, add back a factor on g as a second dial |
| Scrap capacity counts one panel per scrap and fits the panel's footprint ellipse rather than the hexagon, so it is slightly conservative | Use it only for sizing; the exact fit test decides the assignment |
| Bow-tie panels in saddle regions fit scraps poorly | Report fit rate per curvature region; consider a smaller η there |
| The fit test is slow (hundreds of panels × scraps × rotations) | Filter by area and inradius; cache per panel shape |

## Stretch goals

- **Fit-aware planarization.** Add a term to [planarize_hex.m](../matlab/planarize_hex.m) that keeps
  each panel inside its assigned scrap, with the panel's placement in the scrap as a variable. This
  builds the Fit objective into the geometry itself. It is the strongest next contribution.
- **Scrap elongation drives panel shape.** A per-face η can already stretch a panel as far as the
  curvature allows, so a face given a long board could take a long panel in parabolic regions
  (3.5–5:1 boards use 53% of their area with a 5:1 panel, 16% with a round one). Flat regions
  would need a separate cap along each principal direction.
- **Comparison with tangent plane intersection**, the proposal's original route.
- **Car and Three holes.** Leave them as a documented limitation unless time is left over.

## Open decisions

- Final due date and presentation date.
- Stock source: synthetic only, or also a real inventory (for example, measured or photographed
  offcuts)?
- Cutting margin: the kerf to add around each panel before the fit test.
- Rotation: are panels free to rotate within a scrap, or does grain direction matter? Default: free.
