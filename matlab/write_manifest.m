function write_manifest(webdata)
%WRITE_MANIFEST Merge web/data/<id>.meta.json files into web/data/manifest.json.
if nargin < 1
    webdata = fullfile(startup_phcpf(), 'web', 'data');
end
cat = model_catalog();
models = {};
for i = 1:numel(cat)
    metaPath = fullfile(webdata, [cat(i).id '.meta.json']);
    if ~exist(metaPath, 'file'), continue; end
    models{end+1} = jsondecode(fileread(metaPath)); %#ok<AGROW>
end
% The paper's planarity colormap (official repo, code/planarity_colormap.mat), 256 RGB rows,
% thinned to 32 stops for the viewer's legend and face colors.
cm = load('planarity_colormap.mat');
stops = round(linspace(1, size(cm.cm1, 1), 32));
solver = 'MATLAB lsqnonlin (Levenberg-Marquardt), as in the official code';
manifest = struct('generated', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), ...
                  'solver', solver, ...
                  'planarityColormap', round(cm.cm1(stops, :), 4), ...
                  'models', {models});
write_json(fullfile(webdata, 'manifest.json'), manifest);
fprintf('manifest: %d models\n', numel(models));
end
