function root = startup_phcpf()
%STARTUP_PHCPF Put the ScrapShell MATLAB code and the vendored PH-CPF code on the path.
%   Returns the ScrapShell repo root.
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
phcpf = fullfile(root, 'third_party', 'PH-CPF');
addpath(fullfile(phcpf, 'code'));
addpath(fullfile(phcpf, 'data'));
addpath(fullfile(phcpf, 'external', 'geom3d'));
addpath(fullfile(phcpf, 'external', 'WOBJ_toolbox_Version2b'));
addpath(here);  % last added = first on path, so our shims (fast_sparse) win
% The official field optimization and integer rounding call lsqnonlin / lsqlin.
if exist('lsqnonlin', 'file') == 0
    warning('ScrapShell:toolbox', ['The Optimization Toolbox is required (lsqnonlin). ' ...
        'Install it from Home > Add-Ons > Get Add-Ons.']);
end
end
