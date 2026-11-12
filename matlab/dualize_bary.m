function [Vd, Fd, info] = dualize_bary(V, F, strip)
%DUALIZE_BARY Barycentric dual of a polygon mesh (replaces PH-CPF's Dualizer_bin).
%   [Vd, Fd, info] = dualize_bary(V, F, strip)
%   V: nv x 3 primal vertices; F: cell array of 1-based primal faces (consistently oriented).
%   Every interior primal vertex becomes a dual face through the centroids of its incident primal
%   faces, ordered around the vertex (paper: p_bar = (p + p_i + p_i+1)/3 for triangles).
%   strip (default 1): also drop dual faces whose one-ring touches the boundary, like the official
%   dualizer's CLIPPED mode; the primal faces clipped by the mesher near the boundary make those
%   panels irregular.
%   info.primalVertex(k) is the primal vertex of dual face k; info.primalFace(j) the primal face
%   whose centroid is dual vertex j.
if nargin < 3, strip = 1; end
nv = size(V, 1);
nf = numel(F);

% directed half-edge -> face, and one incident face per vertex
deg = cellfun(@numel, F);
I = zeros(sum(deg), 1); J = I; K = I; o = 0;
vf = zeros(nv, 1);
for f = 1:nf
    p = F{f}; d = numel(p);
    I(o+1:o+d) = p; J(o+1:o+d) = p([2:d 1]); K(o+1:o+d) = f;
    vf(p) = f;
    o = o + d;
end
H = sparse(I, J, K, nv, nv);
if any(nonzeros(H) ~= round(nonzeros(H))) || nnz(H) ~= numel(I)
    error('dualize_bary: mesh has repeated directed edges (non-manifold or inconsistently oriented)');
end
% boundary vertices: directed edge without its twin
Hb = (H ~= 0) & (H' == 0);
[bi, bj] = find(Hb);
isB = false(nv, 1); isB([bi; bj]) = true;

% face centroids
C = zeros(nf, 3);
for f = 1:nf, C(f, :) = mean(V(F{f}, :), 1); end

% walk the face ring of each interior vertex
ring = cell(nv, 1);
for v = 1:nv
    if isB(v) || vf(v) == 0, continue; end
    f0 = vf(v); f = f0; r = zeros(1, 0);
    for guard = 1:64
        r(end+1) = f; %#ok<AGROW>
        p = F{f}; d = numel(p); j = find(p == v, 1);
        prev = p(mod(j - 2, d) + 1);
        g = H(v, prev);       % face across edge prev->v contains v->prev
        if g == 0, r = []; break; end
        f = g;
        if f == f0, break; end
    end
    ring{v} = r;
end

% optional stripping of the boundary layer
keep = ~isB & ~cellfun(@isempty, ring);
for s = 1:strip
    near = false(nv, 1);
    near([bi; bj]) = true;
    if s > 1, near = near | ~keep; end
    A = sparse(I, J, 1, nv, nv); A = (A + A') > 0;
    touch = (A * double(near)) > 0;
    keep = keep & ~touch;
    if s == 1, firstLayer = touch; end %#ok<NASGU>
end

% assemble dual, orienting each face like the primal surface around its vertex
pv = find(keep);
used = unique([ring{pv}]);
map = zeros(nf, 1); map(used) = 1:numel(used);
Vd = C(used, :);
Fd = cell(numel(pv), 1);
for k = 1:numel(pv)
    v = pv(k); r = ring{v};
    n = zeros(1, 3);
    for f = r
        p = V(F{f}, :);
        n = n + cross(p(2, :) - p(1, :), p(3, :) - p(1, :));
    end
    q = C(r, :); c = mean(q, 1); m = zeros(1, 3);
    for j = 1:numel(r)
        m = m + cross(q(j, :) - c, q(mod(j, numel(r)) + 1, :) - c);
    end
    if dot(m, n) < 0, r = fliplr(r); end
    Fd{k} = map(r)';
end
info = struct('primalVertex', pv, 'primalFace', used(:), 'boundaryPrimal', isB);
end
