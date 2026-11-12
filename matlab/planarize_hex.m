function [V, hist] = planarize_hex(V0, F, X, T, varargin)
%PLANARIZE_HEX Planarize a polygon mesh near the input surface (replaces PH-CPF's Planarize_bin).
%   [V, hist] = planarize_hex(V0, F, X, T, 'mp', 1)
%   Follows the supplement, Eq. 23: unknowns are vertex positions v and face normals n_f;
%     lambda_p * sum_f sum_{e in f} <n_f, e>^2        planarity (penalty, weight grows per pass)
%     lambda_d * sum_i |v_i - vhat_i|^2               distance to the projection of v~_i on the input
%     lambda_s * sum |(c_f - v_i) - (v_{i+|f|/2} - c_f)|^2   central symmetry (even-degree faces)
%     lambda_c * sum_i <n_v~, v_i - v~_i>^2           normal drift w.r.t. the previous pass
%     lambda_l * sum phi(|e|) + phi(|diagonal|)       barrier against collapsing edges / diagonals
%   with n_f entering the planarity term only through its direction n_f/|n_f| (the supplement's
%   homogeneous parametrization) and |n_f|^2 = 1 as a gauge term. phi is the barrier of paper Eq. 14
%   (LIM): phi(x) = 1/g(x/s) - 1 with
%   g(t) = t^3 - 3t^2 + 3t for 0 < x < s, 0 for x >= s, with s = 0.1 x the length at the start of
%   the pass (as the paper's s = 0.1 s_prev). Each outer pass runs Levenberg-Marquardt (lsqnonlin,
%   250 iterations like the supplement's Ceres solve), then updates v~, c_f, n_v~, vhat, doubles
%   lambda_p and lambda_l; stops once max planarity (Eq. 29) <= mp percent or after 100 passes
%   (Planarize_bin's default). Deviation: the supplement doubles lambda_l every 5 solver
%   iterations; LM cannot change weights mid-solve, so it is doubled per pass.
%   Weights default to the paper's: lambda_p = lambda_c = 0.01, lambda_l = 0.001, lambda_d = 0.5,
%   lambda_s = 1. Coordinates are scaled to unit mean edge length while optimizing.
%   hist.V{k}, hist.planarity{k}, hist.pass(k): snapshots after each pass (k = 1 is the input).
ip = inputParser;
ip.addParameter('mp', 1);
ip.addParameter('maxOuter', 100);
ip.addParameter('maxInner', 250);
ip.addParameter('lp', 0.01); ip.addParameter('lc', 0.01); ip.addParameter('ll', 0.001);
ip.addParameter('ld', 0.5);  ip.addParameter('ls', 1);    ip.addParameter('ln', 1);
ip.addParameter('lpGrowth', 2);
ip.addParameter('lpMax', 1e6);        % numerical guard; the supplement does not say how far lambda_p grows
ip.addParameter('barrierFrac', 0.1);
ip.addParameter('verbose', true);
ip.parse(varargin{:});
o = ip.Results;

n = size(V0, 1); m = numel(F);
deg = cellfun(@numel, F);

% --- normalize units: unit mean edge length, centered
E = edge_list(F);
cen = mean(V0, 1);
s = 1 / mean(sqrt(sum((V0(E(:, 1), :) - V0(E(:, 2), :)).^2, 2)));
Vs = (V0 - cen) * s;  Xs = (X - cen) * s;

% --- static index structures
% face edges for planarity: (face, a, b)
FE = zeros(sum(deg), 3); o2 = 0;
for f = 1:m
    p = F{f}; d = numel(p);
    FE(o2+1:o2+d, :) = [repmat(f, d, 1), p(:), p([2:d 1])'];
    o2 = o2 + d;
end
% symmetry pairs for even faces: (face, i, opposite)
SP = zeros(0, 3);
for f = 1:m
    p = F{f}; d = numel(p);
    if mod(d, 2) == 0
        h = d / 2;
        SP = [SP; repmat(f, d, 1), p(:), p(mod((0:d-1) + h, d) + 1)']; %#ok<AGROW>
    end
end
% diagonals (opposite vertices) for the barrier
DG = unique(sort(SP(:, 2:3), 2), 'rows');
LB = [E; DG];                                     % segments with a length barrier
len0 = sqrt(sum((Vs(LB(:, 1), :) - Vs(LB(:, 2), :)).^2, 2));
if any(len0 < 1e-9)
    warning('planarize_hex:degenerate', '%d zero-length edges/diagonals in the input panels', sum(len0 < 1e-9));
end
% vertex -> faces (for n_v~)
VFi = cell2mat(arrayfun(@(f) F{f}(:), (1:m)', 'UniformOutput', false));
VFf = repelem((1:m)', deg);
VF = sparse(VFi, VFf, 1, n, m);

% --- initial normals: best-fit plane per face, oriented with the winding
N = zeros(m, 3);
for f = 1:m
    q = Vs(F{f}, :); c = mean(q, 1);
    [~, ~, W] = svd(q - c, 'econ'); nf = W(:, 3)';
    w = [0 0 0];
    for j = 1:size(q, 1), w = w + cross(q(j, :) - c, q(mod(j, size(q, 1)) + 1, :) - c); end
    if dot(nf, w) < 0, nf = -nf; end
    N(f, :) = nf;
end

hist = struct('V', {{V0}}, 'planarity', {{planarity_eq29(V0, F)}}, 'pass', 0, 'lp', o.lp);
lp = o.lp; ll = o.ll;
opts = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'SpecifyObjectiveGradient', true, ...
    'Display', 'off', 'MaxIterations', o.maxInner, 'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);
Vt = Vs;
for pass = 1:o.maxOuter
    % per-pass constants: previous positions, barycenters, vertex normals, projections
    C = cell2mat(cellfun(@(f) mean(Vt(f, :), 1), F, 'UniformOutput', false));
    Nv = VF * N; Nv = Nv ./ max(sqrt(sum(Nv.^2, 2)), 1e-12);
    Vhat = closest_point_mesh(Vt, Xs, T);
    sBar = o.barrierFrac * sqrt(sum((Vt(LB(:, 1), :) - Vt(LB(:, 2), :)).^2, 2));   % s = 0.1 s_prev
    w = struct('p', sqrt(lp), 'n', sqrt(o.ln), 'd', sqrt(o.ld), 's', sqrt(o.ls), 'c', sqrt(o.lc), 'l', sqrt(ll));
    fun = @(x) residuals(x, n, m, FE, SP, LB, sBar, C, Nv, Vt, Vhat, w);
    x0 = [reshape(Vt', [], 1); reshape(N', [], 1)];
    x = lsqnonlin(fun, x0, [], [], opts);
    Vt = reshape(x(1:3*n), 3, [])';
    N = reshape(x(3*n+1:end), 3, [])';
    N = renormalize(N, Vt, F);
    Vout = Vt / s + cen;
    pl = planarity_eq29(Vout, F);
    hist.V{end+1} = Vout; hist.planarity{end+1} = pl; hist.pass(end+1) = pass; hist.lp(end+1) = lp;
    if o.verbose
        fprintf('  planarize pass %2d: lambda_p %.3g  max %.3f%%  mean %.3f%%\n', pass, lp, max(pl), mean(pl));
    end
    if max(pl) <= o.mp, break; end
    lp = min(lp * o.lpGrowth, o.lpMax);
    ll = ll * 2;
end
V = Vt / s + cen;
end

function N = renormalize(N, V, F)
% Unit normals for the next pass; a (near-)zero normal is replaced by the panel's best-fit plane
% normal, oriented like the panel's winding.
len = sqrt(sum(N.^2, 2));
bad = find(len < 1e-8);
N = N ./ max(len, 1e-12);
for f = bad'
    q = V(F{f}, :); c = mean(q, 1);
    [~, ~, W] = svd(q - c, 'econ'); nf = W(:, 3)';
    w = [0 0 0];
    for j = 1:size(q, 1), w = w + cross(q(j, :) - c, q(mod(j, size(q, 1)) + 1, :) - c); end
    if dot(nf, w) < 0, nf = -nf; end
    N(f, :) = nf;
end
end

function [r, J] = residuals(x, n, m, FE, SP, LB, sBar, C, Nv, Vt, Vhat, w)
V = reshape(x(1:3*n), 3, [])';
N = reshape(x(3*n+1:end), 3, [])';
vid = @(i) 3 * (i - 1);                % offset of vertex i in x
nid = @(f) 3 * n + 3 * (f - 1);        % offset of normal f in x
R = {}; I = {}; Jc = {}; S = {}; row = 0;

% planarity: w.p * <n_f/|n_f|, v_b - v_a>. Only the direction of n_f matters (the supplement's
% "homogeneous parametrization" of the unit normal). With <n_f, e> and |n_f| = 1 as a mere penalty,
% a growing lambda_p makes collapsing n_f to 0 cheaper than flattening a hard panel; that panel then
% loses all planarity force (and 0/0 at renormalization produced NaN starts).
f = FE(:, 1); a = FE(:, 2); b = FE(:, 3); k = numel(f);
e = V(b, :) - V(a, :); nf = N(f, :);
nn = max(sqrt(sum(nf.^2, 2)), 1e-12);
nh = nf ./ nn;
ne = sum(nh .* e, 2);
R{end+1} = w.p * ne;
dn = w.p * (e - ne .* nh) ./ nn;                 % d/dn_f of <n_f/|n_f|, e>
rr = row + (1:k)';
for c = 1:3
    I{end+1} = rr; Jc{end+1} = nid(f) + c; S{end+1} = dn(:, c);              %#ok<AGROW>
    I{end+1} = rr; Jc{end+1} = vid(b) + c; S{end+1} = w.p * nh(:, c);       %#ok<AGROW>
    I{end+1} = rr; Jc{end+1} = vid(a) + c; S{end+1} = -w.p * nh(:, c);      %#ok<AGROW>
end
row = row + k;

% unit normals: w.n * (|n_f|^2 - 1)
R{end+1} = w.n * (sum(N.^2, 2) - 1);
rr = row + (1:m)';
for c = 1:3
    I{end+1} = rr; Jc{end+1} = nid((1:m)') + c; S{end+1} = 2 * w.n * N(:, c); %#ok<AGROW>
end
row = row + m;

% distance to surface: w.d * (v - vhat), 3 rows per vertex
D = w.d * (V - Vhat);
R{end+1} = reshape(D', [], 1);
rr = row + (1:3*n)';
I{end+1} = rr; Jc{end+1} = (1:3*n)'; S{end+1} = w.d * ones(3*n, 1);
row = row + 3 * n;

% symmetry: w.s * (2 c_f - v_i - v_opp), 3 rows per pair (c_f fixed within the pass)
if ~isempty(SP)
    f = SP(:, 1); i = SP(:, 2); j = SP(:, 3); k = numel(f);
    Rs = w.s * (2 * C(f, :) - V(i, :) - V(j, :));
    R{end+1} = reshape(Rs', [], 1);
    for c = 1:3
        rr = row + 3 * (0:k-1)' + c;
        I{end+1} = rr; Jc{end+1} = vid(i) + c; S{end+1} = -w.s * ones(k, 1); %#ok<AGROW>
        I{end+1} = rr; Jc{end+1} = vid(j) + c; S{end+1} = -w.s * ones(k, 1); %#ok<AGROW>
    end
    row = row + 3 * k;
end

% normal drift: w.c * <n_v~, v - v~>
R{end+1} = w.c * sum(Nv .* (V - Vt), 2);
rr = row + (1:n)';
for c = 1:3
    I{end+1} = rr; Jc{end+1} = vid((1:n)') + c; S{end+1} = w.c * Nv(:, c); %#ok<AGROW>
end
row = row + n;

% length barrier (paper Eq. 14): energy w.l^2 * phi(|seg|), residual w.l * sqrt(phi),
% phi = 1/g - 1, g(t) = t^3 - 3t^2 + 3t with t = |seg|/s, active for |seg| < s.
% Near t = 1, phi ~ (1-t)^3, so sqrt(phi) is C1 and the residual switches on smoothly.
a = LB(:, 1); b = LB(:, 2); k = numel(a);
dv = V(a, :) - V(b, :); L = sqrt(sum(dv.^2, 2));
act = L < sBar & sBar > 0;
r = zeros(k, 1); dr = zeros(k, 1);
t = max(L(act) ./ sBar(act), 1e-9);              % guard: a collapsed segment gives a huge, finite value
g = t.^3 - 3 * t.^2 + 3 * t;
phi = 1 ./ g - 1;
dphi = -(3 * t.^2 - 6 * t + 3) ./ g.^2 ./ sBar(act);   % d phi / d|seg|
sq = sqrt(max(phi, 0));
r(act) = sq;
dr(act) = dphi ./ (2 * max(sq, 1e-150)) .* (sq > 0);
R{end+1} = w.l * r;
g = (w.l * dr ./ max(L, 1e-12)) .* dv;           % d/dv_a; d/dv_b = -g
rr = row + (1:k)';
for c = 1:3
    I{end+1} = rr; Jc{end+1} = vid(a) + c; S{end+1} = g(:, c);  %#ok<AGROW>
    I{end+1} = rr; Jc{end+1} = vid(b) + c; S{end+1} = -g(:, c); %#ok<AGROW>
end
row = row + k;

r = vertcat(R{:});
if nargout > 1
    J = sparse(vertcat(I{:}), vertcat(Jc{:}), vertcat(S{:}), row, numel(x));
end
end

function E = edge_list(F)
E = zeros(0, 2);
for f = 1:numel(F)
    p = F{f}; E = [E; p(:), p([2:end 1])']; %#ok<AGROW>
end
E = unique(sort(E, 2), 'rows');
end
