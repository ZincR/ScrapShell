function out = run_phcpf(id, varargin)
%RUN_PHCPF Run the PH-CPF baseline on one model and export every stage for the web viewer.
%   out = run_phcpf('metro')                 % all stages
%   out = run_phcpf('metro', 'upto', 4)      % stop after a stage (1 = input ... 8 = planarized)
%
%   Writes web/data/<id>.json (+ <id>.meta.json, manifest.json) and results/<id>/ (intermediates).
%   Stages: 1 input, 2 curvature regions, 3 guiding field + sizing ellipses, 4 CPF field (LM
%   iterates), 5 seamless parameterization, 6 primal triangle mesh, 7 hex dual, 8 planarized.
%   Parameters default to runme.m of the official code (n = 6, delta = 0.5%, eta = 0.2%, mp = 1%).
root = startup_phcpf();
ip = inputParser;
ip.addParameter('upto', 8);
ip.addParameter('delta', 0.005);   % max distance to surface, fraction of bbox diagonal
ip.addParameter('eta', 0.002);     % max element area, fraction of total area
ip.addParameter('mp', 1);          % planarity target, percent
ip.addParameter('iter_opt', 0);    % continuity-weight continuation passes (runme retries 1..3)
ip.addParameter('strip', 1);       % boundary layers removed by the dualizer (open shells)
ip.addParameter('maxFrames', 36);  % animation frames exported per animated stage
ip.parse(varargin{:});
opt = ip.Results;

cat = model_catalog();
entry = cat(strcmp({cat.id}, id));
if isempty(entry), error('run_phcpf: unknown model %s', id); end

resdir = fullfile(root, 'results', id);
webdata = fullfile(root, 'web', 'data');
if ~exist(resdir, 'dir'), mkdir(resdir); end
offPath = fullfile(root, 'third_party', 'PH-CPF', 'data', [id '.off']);
mesherBin = fullfile(root, 'build', 'phcpf_mesher', 'phcpf_mesher');

out = struct('id', id, 'label', entry.label, 'group', entry.group, ...
             'params', struct('n', 6, 'delta', opt.delta, 'eta', opt.eta, 'mp', opt.mp, 'iter_opt', opt.iter_opt), ...
             'stages', struct());
meta = struct('id', id, 'label', entry.label, 'group', entry.group, 'file', [id '.json'], ...
              'stages', [], 'metrics', struct(), 'error', '');

try
    %% Stage 1: input surface
    [X, T] = MESH_IO.roff(offPath);
    bbd = norm(max(X) - min(X));
    out.stages.input = struct('V', json_array(X), 'F', json_array(T - 1, 'i32'));
    out.bbd = bbd;
    meta.nv = size(X, 1); meta.nf = size(T, 1); meta.bbd = bbd;
    done(1);
    fprintf('[%s] 1 input: %d vertices, %d faces\n', id, size(X, 1), size(T, 1));
    if opt.upto < 2, finish(); return; end

    %% Stages 2-3: curvature, guiding field, sizing (CPF.set_param, as runme.m)
    params = {'n', 6, 'delta', opt.delta, 'eta', opt.eta, 'mp', opt.mp, ...
              'lsm', 1, 'lc', 1, 'lsi', 1, 'la', 1, 'lin', 1, 'lcj', 1, 'lo', 1, ...
              'siz', 'k', 'optinit', 'k', 'iter_opt', opt.iter_opt, ...
              'datadir', [fileparts(offPath) filesep], 'outdir', [resdir filesep]};
    t0 = tic;
    cpf = CPF(id);
    cpf.set_param(params);
    meta.metrics.timeSetup = toc(t0);
    M = cpf.M; nf = M.nf;

    % regions, with the thresholds of SHAPEOP.guiding_field2 (rho_min 0.01, delta_phip 0.05)
    [~, ~, ~, damin, ~, RHO, PHI] = SHAPEOP.guiding_field2(M);
    region = ones(nf, 1);                                 % 1 elliptic (incl. umbilic)
    region(abs(PHI - pi/4) < 0.05) = 2;                   % 2 parabolic
    region(PHI < pi/4 - 0.05) = 3;                        % 3 hyperbolic
    region(RHO <= 0.01) = 0;                              % 0 planar
    out.stages.curvature = struct('region', json_array(region, 'u8'), 'dmin', json_array(damin));
    meta.metrics.regions = struct('planar', sum(region == 0), 'elliptic', sum(region == 1), ...
                                  'parabolic', sum(region == 2), 'hyperbolic', sum(region == 3));
    done(2);
    fprintf('[%s] 2 curvature: %s\n', id, jsonencode(meta.metrics.regions));
    if opt.upto < 3, finish(); return; end

    % guiding 2-RoSy field (optimized_guiding_field) and sizing ellipses {w : w' g w = 1}
    w = cpf.w; w = w ./ max(sqrt(sum(w.^2, 2)), 1e-12);
    [~, S] = SHAPEOP.curvature_sizing(M, cpf.k_eps, cpf.delta);
    a1 = S.dminf .* sqrt(cpf.delta ./ S.kminfs);
    a2 = S.dmaxf .* sqrt(cpf.delta ./ S.kmaxfs);
    out.stages.guide = struct('dir', json_array(w), 'axis1', json_array(a1), 'axis2', json_array(a2));
    done(3);
    fprintf('[%s] 3 guiding field + sizing ellipses\n', id);
    if opt.upto < 4, finish(); return; end

    %% Stage 4: CPF optimization (Levenberg-Marquardt), recording iterates
    t0 = tic;
    phcpf_recorder('start');
    cleanRec = onCleanup(@() phcpf_recorder('stop'));
    cpf.opt_power_cfs();
    phcpf_recorder('stop');
    meta.metrics.timeField = toc(t0);
    [XS, iters, resn] = phcpf_recorder('get');
    nie = M.nie;
    pick = frame_pick(numel(XS), opt.maxFrames);
    frames = cell(numel(pick), 1);
    for k = 1:numel(pick)
        x = double(XS{pick(k)});
        Uk = reshape(M.EBI * x(2*nie+1 : 2*nie+2*nf), [], 3);
        Vk = reshape(M.EBI * x(2*nie+2*nf+1 : end), [], 3);
        frames{k} = struct('U', json_array(Uk), 'V', json_array(Vk));
    end
    [~, Ssing] = cpf.compute_matching_and_sing(cpf.U, cpf.V);
    sv = find(Ssing); si = Ssing(sv);
    lico = dot(cross(cpf.U, cpf.V, 2), M.Nf, 2);
    out.stages.field = struct('frames', {frames}, 'iters', iters(pick), 'energy', resn(pick), ...
        'sing', struct('v', json_array(sv - 1, 'i32'), 'idx', json_array(si, 'i32')));
    meta.metrics.nSing = numel(sv);
    meta.metrics.licoMin = min(lico) / mean(abs(lico));
    meta.metrics.lmIterations = numel(XS) - 1;
    done(4);
    fprintf('[%s] 4 CPF field: %d LM iterates, %d singularities, LICO min/mean %.3f, %.1fs\n', ...
        id, numel(XS), numel(sv), meta.metrics.licoMin, meta.metrics.timeField);
    clear XS
    if opt.upto < 5, finish(); return; end

    %% Stages 5-6: seamless integration + meshing (tools/phcpf_mesher, Directional)
    [gu, gv] = CPF.uv2grads(M.EB * cpf.U(:), M.EB * cpf.V(:), nf);
    gradu = reshape(M.EBI * gu, [], 3); gradv = reshape(M.EBI * gv, [], 3);
    Fraw = cpf.grads2rawfield(gradu, gradv);
    rawPath = fullfile(resdir, [id '.rawfield']);
    MESH_IO.wrawfield(rawPath, Fraw);
    prefix = fullfile(resdir, id);
    % As in CPF.uv2mesh: MATLAB's matching and singularities -> cut + comb + seamless Poisson system
    % (phcpf_mesher setup = MeshCutter + Parametrization) -> official BarrierPoissonSolve.m rounding
    % -> isoline meshing of the rounded functions (phcpf_mesher mesh = Mesher).
    write_cpf_matching(cpf, prefix);
    inArgs = sprintf('--matching "%s_matching.txt" --sing "%s_singin.txt"', prefix, prefix);
    t0 = tic;
    run_mesher(mesherBin, 'setup', offPath, rawPath, prefix, inArgs);
    meta.metrics.timeSetupParam = toc(t0);
    t0 = tic;
    [xr, ok, binfo] = barrier_rounding(cpf, prefix);
    meta.metrics.timeRounding = toc(t0);
    meta.metrics.poissonError = binfo.poissonError;
    meta.metrics.nIntegerVars = binfo.nInteger;
    meta.metrics.flippedAfterRounding = binfo.flippedFaces;
    fprintf('[%s] 5 BarrierPoissonSolve: %d integer variables, Poisson error %.2g, success %d, %d flipped faces, %.1fs\n', ...
        id, binfo.nInteger, binfo.poissonError, ok, binfo.flippedFaces, meta.metrics.timeRounding);
    if ok ~= 1
        error('BarrierPoissonSolve failed, cannot find integer param.');   % CPF.param_and_mesh message
    end
    writematrix(xr(:), [prefix '_x.txt'], 'FileType', 'text');
    t0 = tic;
    run_mesher(mesherBin, 'mesh', offPath, rawPath, prefix, sprintf('%s --x "%s_x.txt"', inArgs, prefix));
    meta.metrics.timeMesher = toc(t0);
    C = read_corner([prefix '_corner.txt']);           % nf x 18 (3 corners x 6 functions)
    H = [C(:, 1:3), C(:, 7:9), C(:, 13:15)];          % h1..h3 at each corner
    SM = readmatrix([prefix '_sing.txt'], 'FileType', 'text');   % singularities used for cutting
    if isempty(SM), SM = zeros(0, 2); end
    out.stages.param = struct('h', json_array(H), ...
        'sing', struct('v', json_array(SM(:, 1), 'i32'), 'idx', json_array(SM(:, 2), 'i32')));
    ms = jsondecode(fileread([prefix '_stats.json']));
    meta.metrics.nSingMesher = ms.nSing;
    if ms.nF == 0
        error('phcpf_mesher produced no polygons from the rounded parameterization');
    end
    done(5);
    [Vp, Fp] = read_poly_off([prefix '_tri.off']);
    [Dp, Ip] = flatten_faces(Fp);
    out.stages.primal = struct('V', json_array(Vp), 'D', json_array(Dp, 'i32'), 'idx', json_array(Ip - 1, 'i32'));
    meta.metrics.nPrimal = numel(Fp);
    meta.metrics.primalDegrees = accumarray(Dp, 1)';
    done(6);
    fprintf('[%s] 5-6 param + primal mesh: %d polygons (degrees %s), mesher singularities %d vs field %d\n', ...
        id, numel(Fp), mat2str(meta.metrics.primalDegrees), ms.nSing, meta.metrics.nSing);
    if opt.upto < 7, finish(); return; end

    %% Stage 7: barycentric hex dual
    [Vd, Fd] = dualize_bary(Vp, Fp, opt.strip);
    [Dd, Id] = flatten_faces(Fd);
    pd = planarity_eq29(Vd, Fd);
    write_poly_off(fullfile(resdir, [id '_dual.off']), Vd, Fd);
    out.stages.dual = struct('V', json_array(Vd), 'D', json_array(Dd, 'i32'), 'idx', json_array(Id - 1, 'i32'), ...
        'planarity', json_array(pd));
    meta.metrics.planMaxBefore = max(pd); meta.metrics.planMeanBefore = mean(pd);
    done(7);
    fprintf('[%s] 7 dual: %d panels, planarity max %.2f%% mean %.3f%%\n', id, numel(Fd), max(pd), mean(pd));
    if opt.upto < 8, finish(); return; end

    %% Stage 8: planarization
    t0 = tic;
    [Vq, hist] = planarize_hex(Vd, Fd, X, T, 'mp', opt.mp);
    meta.metrics.timePlanar = toc(t0);
    write_poly_off(fullfile(resdir, [id '_planar.off']), Vq, Fd);
    pick = frame_pick(numel(hist.V), opt.maxFrames);
    pf = cell(numel(pick), 1);
    for k = 1:numel(pick)
        pf{k} = struct('V', json_array(hist.V{pick(k)}), 'planarity', json_array(hist.planarity{pick(k)}));
    end
    out.stages.planar = struct('frames', {pf}, 'iters', hist.pass(pick));
    km = metrics_phcpf(Vq, Fd, X, T);
    for fn = fieldnames(km)', meta.metrics.(fn{1}) = km.(fn{1}); end
    meta.metrics.planarPasses = hist.pass(end);
    done(8);
    fprintf('[%s] 8 planarized: max %.3f%% mean %.3f%%, dist max %.4f bbd, %d passes, %.1fs\n', ...
        id, km.planMax, km.planMean, km.distMax, hist.pass(end), meta.metrics.timePlanar);
catch err
    meta.error = err.message;
    fprintf(2, '[%s] stopped: %s\n', id, getReport(err, 'basic'));
end
finish();

    function done(k)
        meta.stages(end+1) = k;
    end
    function finish()
        write_json(fullfile(webdata, [id '.json']), out);
        write_json(fullfile(webdata, [id '.meta.json']), meta);
        write_manifest(webdata);
    end
end

function run_mesher(bin, mode, offPath, rawPath, prefix, extra)
% Run tools/phcpf_mesher. The isoline mesher can loop forever on a degenerate parameterization, so
% each call is capped at 120 s (perl alarm; macOS has no `timeout`) and the attempt then fails.
cmd = sprintf(['env -u DYLD_LIBRARY_PATH -u DYLD_FALLBACK_LIBRARY_PATH ' ...
    'perl -e ''alarm shift @ARGV; exec @ARGV'' 120 "%s" %s "%s" "%s" "%s" %s'], ...
    bin, mode, offPath, rawPath, prefix, extra);
[status, msg] = system(cmd);
fprintf('%s', msg);
if status ~= 0, error('phcpf_mesher %s failed (%d): %s', mode, status, msg); end
end

function pick = frame_pick(n, maxFrames)
% Frames to export: all if few, otherwise denser early (the field changes most at the start).
if n <= maxFrames, pick = 1:n; return; end
pick = unique(round(1 + (n - 1) * linspace(0, 1, maxFrames).^2));
end

function C = read_corner(path)
fid = fopen(path, 'r');
hdr = fscanf(fid, '%d', 2);
C = fscanf(fid, '%f', [3 * hdr(2), hdr(1)])';
fclose(fid);
end

function [D, idx] = flatten_faces(F)
D = cellfun(@numel, F);
idx = cell2mat(cellfun(@(f) f(:), F, 'UniformOutput', false));
end
