# Project plan: Scrap Shells

Proposed 2026-11-21. The final due date is not yet confirmed (see [Open decisions](#open-decisions)).
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

So the extension needs only one new input to the field: a per-face factor $\sigma \le 1$ computed
from the stock. Replacing g by $g/\sigma^2$ on a face shrinks the lattice there by σ, because U and V
map the lattice axes onto the surface. Panels can then be smaller than their footprint, never
larger.

### Pipeline

| Proposal step | Planned implementation | Status |
| --- | --- | --- |
| **Seed:** one footprint per scrap, sized by the stock | Size field σ over the surface, computed by 1D optimal transport between the deposit's scrap sizes and the footprint sizes | To do |
| **Grow:** regions grow toward their scrap's size | PH-CPF field and seamless parameterization with the σ-scaled metric | Baseline done; σ hook to do |
| **Intersect:** flat panels that close the shell | PH-CPF dual and planarization. Closed by construction | Done (baseline) |
| **Match:** pair each scrap with a panel it fits | Exact fit test plus linear assignment between panels and scraps | To do |
| **Check / Revise:** reseed where it fails | Feedback loop: shrink σ where panels fit no scrap | To do |
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
2. **Size field σ by 1D optimal transport.** New MATLAB or Python function.
   - Footprint size per face: the semi-axes are already computed in stage 3 of
     [run_phcpf.m](../matlab/run_phcpf.m) (`a1`, `a2`). Predict the panel's circumradius as c × the
     larger semi-axis, with c calibrated once on a baseline run.
   - Sort faces by footprint and scraps by inradius, both largest first. Walk the two lists together.
     Each scrap claims the surface area of one panel of its size; those faces get
     σ = min(1, scrap inradius ÷ predicted circumradius). This is the monotone (optimal) coupling: big
     scraps go where the surface has room. The deposit holds far more than one surface needs, so
     this is a partial transport: it also chooses which scraps to use.
   - Smooth log σ with a few Laplacian steps and limit its change between neighbouring faces, so the
     field stays integrable.
   - The inradius rule is conservative. In the fit preview, centring a panel on the scrap's incenter
     already passes 1.6–3.4× as many scraps as the inscribed-circle test
     ([stock.md](stock.md#fit-preview)). Calibrate a factor on the scrap side too, or let the
     feedback loop (item 7) recover the slack.
   - *Done when* σ can be plotted on metro.
3. **Hook σ into CPF.**
   - Add a `sigma` option to `CPF.set_param` that applies $g \leftarrow D g D$, with
     $D = \mathrm{diag}(1/\sigma, 1/\sigma)$, right after g is built
     ([CPF.m:302-313](../third_party/PH-CPF/code/CPF.m#L302-L313)).
   - It must go inside `set_param`, not after it: `set_param` caches the mean edge length under g
     (`mel2d`, line 320) and uses it to normalize the continuity and alignment weights (lines 356 and
     362). `opt_init` normalizes the initial field by g, so it follows automatically.
   - Raise η so the stock, not η, caps panel size in flat regions.
   - Pass σ through [run_phcpf.m](../matlab/run_phcpf.m). Record the edit in
     [PATCHES.md](../third_party/PH-CPF/PATCHES.md).
4. **Sanity checks on metro.**
   - Constant σ = 0.7: about twice as many panels, same planarity and distance.
   - A smooth gradient in σ: rounding still succeeds, and the singularity count stays reasonable.
   - *Done when* both runs reach stage 8 and appear in the viewer.

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
7. **Feedback loop.** Shrink σ around panels that fit no scrap (or drop already-used scraps from the
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
| 1 | Nov 21 – Nov 27 (Thanksgiving Nov 26) | Phase 1: stock generator, size field, σ hook, sanity checks |
| 2 | Nov 28 – Dec 4 | Phase 2: fit test, assignment, feedback loop, viewer stage 9 |
| 3 | Dec 5 – Dec 11 | Phase 3 experiments; Phase 4 report and presentation |

The work splits cleanly between two people: Python (stock, fit test, assignment: items 1 and 5–7)
and MATLAB (size field, σ hook, feedback driver: items 2–4).

## Risks

| Risk | Mitigation |
| --- | --- |
| Varying panel size makes the field harder to integrate: more singularities, rounding failures (the same failure as Car and Three holes) | Smooth σ and limit its gradient; reuse the continuation in [run_with_retries.m](../matlab/run_with_retries.m); run experiments on the models that already succeed |
| The grid only roughly hits the target size, so first-pass panels miss their scraps | The feedback loop (item 7) |
| The inradius proxy underestimates what fits in elongated scraps (measured: centring the panel on the incenter passes 1.6–3.4× as many scraps) | Use the proxy only for sizing; the exact fit test decides the assignment |
| Bow-tie panels in saddle regions fit scraps poorly | Report fit rate per curvature region; consider smaller σ there |
| The fit test is slow (hundreds of panels × scraps × rotations) | Filter by area and inradius; cache per panel shape |

## Stretch goals

- **Fit-aware planarization.** Add a term to [planarize_hex.m](../matlab/planarize_hex.m) that keeps
  each panel inside its assigned scrap, with the panel's placement in the scrap as a variable. This
  builds the Fit objective into the geometry itself. It is the strongest next contribution.
- **Scrap elongation drives panel shape.** Long scraps could produce brick-shaped panels, through the
  metric's anisotropy, in parabolic and flat regions.
- **Comparison with tangent plane intersection**, the proposal's original route.
- **Car and Three holes.** Leave them as a documented limitation unless time is left over.

## Open decisions

- Final due date and presentation date.
- Stock source: synthetic only, or also a real inventory (for example, measured or photographed
  offcuts)?
- Cutting margin: the kerf to add around each panel before the fit test.
- Rotation: are panels free to rotate within a scrap, or does grain direction matter? Default: free.
