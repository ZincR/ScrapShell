function run_with_retries(id, varargin)
%RUN_WITH_RETRIES Run run_phcpf, retrying with iter_opt = 1..10 when the pipeline stops early.
%   Paper Sec. 4.5: the continuity weight is doubled (warm-up passes, iter_opt) until integer
%   meshing succeeds, up to 10 times; this makes the field closer to integrable and removes
%   parameter-space flips. If every attempt fails, the last attempt's data (up to the stage it
%   reached) stays in web/data and the viewer shows how far it got.
maxDoublings = 10;
for iter_opt = 0:maxDoublings
    run_phcpf(id, 'iter_opt', iter_opt, varargin{:});
    meta = jsondecode(fileread(fullfile(startup_phcpf(), 'web', 'data', [id '.meta.json'])));
    if max(meta.stages) == 8, return; end
    if iter_opt < maxDoublings
        fprintf('[%s] reached stage %d (%s); retrying with iter_opt = %d\n', ...
            id, max(meta.stages), meta.error, iter_opt + 1);
    else
        fprintf('[%s] gave up after %d doublings; stopped at stage %d (%s)\n', ...
            id, maxDoublings, max(meta.stages), meta.error);
    end
end
end
