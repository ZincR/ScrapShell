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
   - **Fit, not count, is the constraint.** A quick fit preview finds that only 7–38% of scraps
     can hold a typical baseline panel (Metro 7%: about 300 scraps for 436 panels). We keep the
     0.5 m median on purpose: the size field should have real work to do. Revisit once the exact
     fit test (item 5) exists.
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
   - **From target to η.** On each face, solve for $\epsilon_f$ so the panel has the target area (a
     quadratic in $\epsilon_f^2$), then raise $\epsilon_f$ where needed to keep the aspect at most 2
     (also closed form).
   - Smooth log η with a few Laplacian steps and limit its change between neighbouring faces, so the
     field stays integrable.
   - *Done when* η can be plotted on metro.
3. **Hook η into CPF.** Do this and item 4 before the rest of item 2: they test the main risk
   (see [Risks](#risks)).
   - Let `CPF.set_param` take η as one value per face, each a fraction of the total area like the
     single value. Two edits make the sizing code elementwise: `cpf.delta ./ cpf.eta` for `k_eps`
     ([CPF.m:296](../third_party/PH-CPF/code/CPF.m#L296)), and `e.^2` for `e^2` in
     `SHAPEOP.curvature_sizing` and `check_curvature_sizing`
     ([SHAPEOP.m:267-268](../third_party/PH-CPF/code/SHAPEOP.m#L267-L268), 291–292). A single η
     behaves exactly as before.
   - g is built inside `set_param`, so everything that depends on it follows: the mean edge length
     under g (`mel2d`, line 320), the continuity and alignment weights it normalizes (lines 356 and
     362), and `opt_init`, which normalizes the initial field by g.
   - Pass η through [run_phcpf.m](../matlab/run_phcpf.m), as one value or a file with one per face,
     and store the per-face values in the viewer data: [sizing/footprint.py](../sizing/footprint.py)
     recovers curvatures by inverting the stored ellipses with the run's η. Record the edits in
     [PATCHES.md](../third_party/PH-CPF/PATCHES.md).
4. **Sanity checks on metro.**
   - η given per face, equal to the baseline's everywhere: the same result as the baseline run.
   - A constant 4× η: about half as many panels, same planarity and distance.
   - η varying smoothly from 1× to 16× across the model: rounding still succeeds, and the
     singularity count stays reasonable.
   - *Done when* the runs reach stage 8 and appear in the viewer.

### Phase 2: matching, cutting and feedback

5. **Fit test.** Python, with pyclipper or shapely.
   - Flatten each planarized panel P onto its best-fit plane. At ≤ 1% planarity this is safe.
   - P fits scrap S if, for some rotation, there is a translation that puts P inside S. That
     translation set is S shrunk by P: the complement of (S's complement ⊕ −P), a Minkowski sum.
     Sample about 36 rotations.
   - Filter by area and inradius first, so each panel is tested against only a few scraps.
   - Works for non-convex bow-tie panels too.
   - Optionally grow every panel by the saw kerf before testing.
6. **Assignment.** `scipy.optimize.linear_sum_assignment`. The cost of a pair is the waste (scrap
   area minus panel area), with a very large cost for pairs that do not fit. Outputs: the
   assignment, the panels that fit nowhere, utilization, and each panel's placement in its scrap.
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
| Varying panel size makes the field harder to integrate: more singularities, rounding failures (the same failure as Car and Three holes) | Test a varying η on metro first (items 3–4); smooth log η and limit its gradient; reuse the continuation in [run_with_retries.m](../matlab/run_with_retries.m); run experiments on the models that already succeed |
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
