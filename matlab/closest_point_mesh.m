function [Q, dist, fid] = closest_point_mesh(P, X, T, rings)
%CLOSEST_POINT_MESH Closest points on a triangle mesh (X, T) for query points P.
%   Candidates per query: triangles within RINGS (default 2) vertex rings of the nearest mesh
%   vertex (dsearchn). Exact point-triangle projection (Ericson, Real-Time Collision Detection 5.1.5)
%   among the candidates. Accurate for points near the surface, which is all we query.
if nargin < 4, rings = 2; end
nv = size(X, 1); nt = size(T, 1);
VT = sparse(T(:), repmat((1:nt)', 3, 1), 1, nv, nt);         % vertex-triangle incidence
A = double((VT * VT') > 0);                                   % vertex adjacency incl. self
R = A;
for r = 2:rings, R = double((R * A) > 0); end
near = dsearchn(X, P);
m = size(P, 1);
Q = zeros(m, 3); dist = zeros(m, 1); fid = zeros(m, 1);
RT = double((R * VT) > 0);                                    % vertex -> candidate triangles
for i = 1:m
    cand = find(RT(near(i), :));
    [q, d2, j] = point_triangles(P(i, :), X(T(cand, 1), :), X(T(cand, 2), :), X(T(cand, 3), :));
    Q(i, :) = q; dist(i) = sqrt(d2); fid(i) = cand(j);
end
end

function [q, d2, j] = point_triangles(p, a, b, c)
% closest point from p to each triangle (rows of a, b, c); return the best.
n = size(a, 1);
ab = b - a; ac = c - a; ap = p - a;
d1 = sum(ab .* ap, 2); d2v = sum(ac .* ap, 2);
bp = p - b; d3 = sum(ab .* bp, 2); d4 = sum(ac .* bp, 2);
cp = p - c; d5 = sum(ab .* cp, 2); d6 = sum(ac .* cp, 2);
va = d3 .* d6 - d5 .* d4; vb = d5 .* d2v - d1 .* d6; vc = d1 .* d4 - d3 .* d2v;
Qc = zeros(n, 3); done = false(n, 1);
r = d1 <= 0 & d2v <= 0;                         Qc(r, :) = a(r, :); done = done | r;
r = ~done & d3 >= 0 & d4 <= d3;                  Qc(r, :) = b(r, :); done = done | r;
r = ~done & vc <= 0 & d1 >= 0 & d3 <= 0;
t = d1 ./ (d1 - d3);                             Qc(r, :) = a(r, :) + t(r) .* ab(r, :); done = done | r;
r = ~done & d6 >= 0 & d5 <= d6;                  Qc(r, :) = c(r, :); done = done | r;
r = ~done & vb <= 0 & d2v >= 0 & d6 <= 0;
t = d2v ./ (d2v - d6);                           Qc(r, :) = a(r, :) + t(r) .* ac(r, :); done = done | r;
r = ~done & va <= 0 & (d4 - d3) >= 0 & (d5 - d6) >= 0;
t = (d4 - d3) ./ ((d4 - d3) + (d5 - d6));        Qc(r, :) = b(r, :) + t(r) .* (c(r, :) - b(r, :)); done = done | r;
r = ~done;
den = 1 ./ (va + vb + vc); v = vb .* den; w = vc .* den;
Qc(r, :) = a(r, :) + v(r) .* ab(r, :) + w(r) .* ac(r, :);
dd = sum((Qc - p).^2, 2);
[d2, j] = min(dd);
q = Qc(j, :);
end
