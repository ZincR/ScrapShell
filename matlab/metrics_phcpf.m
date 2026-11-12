function k = metrics_phcpf(V, F, X, T)
%METRICS_PHCPF Panel statistics for a polygon mesh (V, F) approximating the input mesh (X, T).
%   planMax / planMean: paper Eq. 29 planarity in percent (max = epsilon_p, mean = mean epsilon_p).
%   distMax / distMean: one-sided distance from the panel mesh (vertices + face centroids) to the
%   input surface, divided by the input bounding-box diagonal (paper's epsilon_d).
deg = cellfun(@numel, F);
p = planarity_eq29(V, F);
C = cell2mat(cellfun(@(f) mean(V(f, :), 1), F, 'UniformOutput', false));
[~, d] = closest_point_mesh([V; C], X, T);
bbd = norm(max(X) - min(X));
k = struct('nPanels', numel(F), 'nHex', sum(deg == 6), 'nOther', sum(deg ~= 6), ...
           'planMax', max(p), 'planMean', mean(p), ...
           'distMax', max(d) / bbd, 'distMean', mean(d) / bbd);
end
