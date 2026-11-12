# Vendored PH-CPF: provenance and local changes

Source: https://github.com/michaled/PH-CPF at commit `0d119e1` (2024-04-12), by Kacper Pluta,
Michal Edelstein, Amir Vaxman and Mirela Ben-Chen. Licensed CC BY-SA 4.0 (see `LICENSE.txt`);
our modifications below are shared under the same license.

## Removed
- `external/{MeshCutter,Parameterization,Mesher,Dualizer,Planarizer}/*.exe` and their DLLs
  (Windows-only). Replaced by `tools/phcpf_mesher` (Directional) and `matlab/dualize_bary.m`,
  `matlab/planarize_hex.m`.
- `external/fast_sparse/fast_sparse.mexw64` (Windows-only). Replaced by `matlab/fast_sparse.m`,
  which forwards to `sparse`.
- `teaser.jpg`.

## Edited
- `code/OPT.m`, `OPT.lsqnonlin`: when `phcpf_recorder('active')` is true, attach an `OutputFcn`
  that stores every Levenberg-Marquardt iterate, so the viewer can animate the field optimization.
  No effect on the optimization itself.
- `code/OPT.m`, `OPT.lsqnonlin`: commented out `opt.CheckGradients = false`. The option no longer
  exists in `optimoptions` (MATLAB R2026a raises an error); it was disabled anyway.

## Unchanged but driven differently
- `code/BarrierPoissonSolve.m` runs unmodified, called from `matlab/barrier_rounding.m` with the
  variables the official `Parametrization_bin` wrote to `poisson.mat` (now exported by
  `tools/phcpf_mesher setup`). Its integer set also contains the singular-vertex variables, which the
  current Directional keeps separate from the seam translations.

## Known upstream discrepancies (left as-is for fidelity)
- `CPF.injectivity` (CPF.m L1066): `gv2` uses `V2(1,:)` where `U2(1,:)` is expected; only matters
  when the LICO barrier is active.
- Sizing: `set_param` uses gamma^2 = delta (no factor 2, no 2*pi in epsilon), whereas the paper's
  Eq. 25 has gamma^2 = 2*delta and epsilon = 2*pi*delta/eta.
- Paper Eq. 16 (smoothness) and Eq. 27 (V0 = -J U0) differ from the code, which follows the
  supplement (Eq. 19) and V0 = J U0.
- `runme.m` retries after a barrier failure without resetting `cpf.lc`, so the continuity weight
  keeps doubling across retries.
