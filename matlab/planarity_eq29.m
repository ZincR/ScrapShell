function p = planarity_eq29(V, F)
%PLANARITY_EQ29 Per-face planarity error in percent (paper Eq. 29), same as MESHP.planarity_general.
%   For every cyclic 4-tuple (v1..v4) of a face: |<n, v2 - v1>| / mean diagonal length, with n the
%   unit normal of the two diagonals (v3-v1) x (v4-v2). Face value = 100 * RMS over its tuples.
%   V: n x 3; F: cell array of 1-based faces. Vectorized per face degree.
nf = numel(F);
p = zeros(nf, 1);
deg = cellfun(@numel, F);
for d = unique(deg)'
    ids = find(deg == d);
    if d < 4, continue; end           % triangles are planar
    M = vertcat(F{ids});              % |ids| x d
    q = zeros(numel(ids), d);
    for j = 1:d
        i1 = M(:, j); i2 = M(:, mod(j, d) + 1); i3 = M(:, mod(j + 1, d) + 1); i4 = M(:, mod(j + 2, d) + 1);
        d1 = V(i3, :) - V(i1, :); d2 = V(i4, :) - V(i2, :);
        x = cross(d1, d2, 2);
        nx = sqrt(sum(x.^2, 2));
        denom = nx .* (sqrt(sum(d1.^2, 2)) + sqrt(sum(d2.^2, 2))) / 2;
        num = sum(x .* (V(i2, :) - V(i1, :)), 2);
        qq = num ./ denom;
        qq(abs(denom) <= 1e-9) = 0;
        q(:, j) = qq;
    end
    p(ids) = 100 * sqrt(sum(q.^2, 2) / d);
end
end
