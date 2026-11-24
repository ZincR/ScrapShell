function run_eta_checks(names)
%RUN_ETA_CHECKS Sanity checks of the per-face eta hook on Metro (docs/plan.md, Phase 1, item 4).
%   run_eta_checks()
%
%   Variant runs of Metro, written beside the baseline in the viewer's "Sizing checks" group. Three checks:
%     metro_eta_1x    eta given per face, equal to the baseline's 0.2% everywhere. The sizing is
%                     the same, so the result must be the baseline's.
%     metro_eta_4x    4x the baseline's eta everywhere, per face: about half as many panels.
%     metro_eta_ramp  eta rising log-linearly from 1x to 16x along Metro's longest axis, passed as a
%                     file (results/metro_eta_ramp.txt) to test that route too.
%   Two diagnostics, to tell panel size from panel stretch:
%     metro_eta_4x_a2 4x eta, lowered where a footprint would be more stretched than 2:1. Needs
%                     results/metro_eta_4x_a2.txt, from
%                     python -m sizing.eta metro --eta 0.008 --max-aspect 2 --out results/metro_eta_4x_a2.txt
%     metro_eta_4x_delta 4x delta and 4x eta: the baseline's footprints, twice as large and no more
%                     stretched (eps = delta / eta is unchanged).
%   Two smoothed size fields (plan item 2, "Smooth the target"): the panel areas of 4x eta, smoothed
%   with 400 Laplacian steps and converted back to eta per face, to test whether a rough metric is
%   what breaks rounding at 4x. They need eta files from sizing/eta.py:
%     metro_eta_4x_smooth     python -m sizing.eta metro --eta 0.008 --smooth 400 --out results/metro_eta_4x_smooth.txt
%     metro_eta_4x_smooth_a2  python -m sizing.eta metro --eta 0.008 --smooth 400 --max-aspect 2 --out results/metro_eta_4x_smooth_a2.txt
%   Each run goes through run_with_retries, as the baseline did. run_eta_checks(names) runs only
%   the named ones.
if nargin < 1
    names = {'metro_eta_1x', 'metro_eta_4x', 'metro_eta_ramp', 'metro_eta_4x_a2', 'metro_eta_4x_delta', ...
             'metro_eta_4x_smooth', 'metro_eta_4x_smooth_a2'};
end
names = cellstr(names);
root = startup_phcpf();
[X, T] = MESH_IO.roff(fullfile(root, 'third_party', 'PH-CPF', 'data', 'metro.off'));
nf = size(T, 1);
base = 0.002;
group = 'Sizing checks';

if ismember('metro_eta_1x', names)
    run_with_retries('metro', 'eta', base * ones(nf, 1), 'out', 'metro_eta_1x', ...
        'label', 'Metro, η per face = baseline', 'group', group);
end
if ismember('metro_eta_4x', names)
    run_with_retries('metro', 'eta', 4 * base * ones(nf, 1), 'out', 'metro_eta_4x', ...
        'label', 'Metro, η 4×', 'group', group);
end
if ismember('metro_eta_ramp', names)
    C = (X(T(:, 1), :) + X(T(:, 2), :) + X(T(:, 3), :)) / 3;
    [~, ax] = max(max(X) - min(X));
    t = (C(:, ax) - min(X(:, ax))) / (max(X(:, ax)) - min(X(:, ax)));
    etaFile = fullfile(root, 'results', 'metro_eta_ramp.txt');
    writematrix(base * 16 .^ t, etaFile, 'FileType', 'text');
    run_with_retries('metro', 'eta', etaFile, 'out', 'metro_eta_ramp', ...
        'label', 'Metro, η 1× to 16×', 'group', group);
end
if ismember('metro_eta_4x_a2', names)
    run_file(root, group, 'metro_eta_4x_a2', '--eta 0.008 --max-aspect 2', 'Metro, η 4×, aspect ≤ 2');
end
if ismember('metro_eta_4x_delta', names)
    run_with_retries('metro', 'delta', 4 * 0.005, 'eta', 4 * base, 'out', 'metro_eta_4x_delta', ...
        'label', 'Metro, δ and η 4×', 'group', group);
end
if ismember('metro_eta_4x_smooth', names)
    run_file(root, group, 'metro_eta_4x_smooth', '--eta 0.008 --smooth 400', 'Metro, η 4×, smoothed');
end
if ismember('metro_eta_4x_smooth_a2', names)
    run_file(root, group, 'metro_eta_4x_smooth_a2', '--eta 0.008 --smooth 400 --max-aspect 2', ...
        'Metro, η 4×, smoothed, aspect ≤ 2');
end
end

function run_file(root, group, oid, args, label)
% A run whose eta comes from results/<oid>.txt, written by sizing/eta.py with these arguments.
etaFile = fullfile(root, 'results', [oid '.txt']);
if ~exist(etaFile, 'file')
    error('run_eta_checks: %s is missing; make it with python -m sizing.eta metro %s --out results/%s.txt', ...
        etaFile, args, oid);
end
run_with_retries('metro', 'eta', etaFile, 'out', oid, 'label', label, 'group', group);
end
