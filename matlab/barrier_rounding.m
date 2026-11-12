function [xcurr, success, info] = barrier_rounding(cpf, prefix)
%BARRIER_ROUNDING Integer rounding with the official BarrierPoissonSolve.m (unmodified).
%   Loads the seamless Poisson system exported by `phcpf_mesher setup` into the variables the official
%   Parametrization_bin wrote to poisson.mat (A, b, C, G, x0, integerIndices, x2CornerMat, V, F, FN, N),
%   then runs the official script in this function's workspace. It rounds the integer variables one
%   at a time, re-solving with Levenberg-Marquardt under a barrier that keeps each face's N gradient
%   vectors in their original cyclic order (no flipped triangles).
%   Adaptation: in the current Directional the singular vertices are separate variables, so they are
%   rounded together with the seam translations (an index +-1/6 cone must sit on a lattice point).
A = read_triplets([prefix '_A.txt']);
b = readmatrix([prefix '_b.txt'], 'FileType', 'text');
C = read_triplets([prefix '_C.txt']);
G = read_triplets([prefix '_G.txt']);
x2CornerMat = read_triplets([prefix '_x2c.txt']);
x0 = readmatrix([prefix '_x0.txt'], 'FileType', 'text');
ints = readmatrix([prefix '_int.txt'], 'FileType', 'text');
sings = readmatrix([prefix '_singidx.txt'], 'FileType', 'text');
integerIndices = unique([ints(:); sings(:)]) + 1;
V = readmatrix([prefix '_cutV.txt'], 'FileType', 'text');
F = readmatrix([prefix '_cutF.txt'], 'FileType', 'text') + 1;
FN = readmatrix([prefix '_FN.txt'], 'FileType', 'text');
meta = jsondecode(fileread([prefix '_meta.json']));
N = meta.N; %#ok<NASGU>
poisson_mat_out = [prefix '_poisson.mat'];
save(poisson_mat_out, 'A', 'b', 'C', 'G', 'x0', 'integerIndices', 'x2CornerMat', 'V', 'F', 'FN', 'N');
clear A b C G x0 integerIndices x2CornerMat V F FN N ints sings meta
xcurr = []; success = 0; minBijectivity = []; fraction = 1; %#ok<NASGU> set by the script below

BarrierPoissonSolve;   % official script: defines xcurr, success, minBijectivity, perr inputs (A, b, x0)

info = struct('nInteger', numel(integerIndices), 'poissonError', norm(A*x0 - b) / norm(A*x0), ...
              'minBijectivity', min(minBijectivity), 'flippedFaces', sum(minBijectivity <= 0), ...
              'fraction', fraction);
end
