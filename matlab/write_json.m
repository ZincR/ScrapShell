function write_json(path, s)
%WRITE_JSON Encode a struct with jsonencode and write it to PATH (creating the folder).
folder = fileparts(path);
if ~isempty(folder) && ~exist(folder, 'dir'), mkdir(folder); end
fid = fopen(path, 'w');
if fid < 0, error('write_json: cannot open %s', path); end
cleaner = onCleanup(@() fclose(fid));
fwrite(fid, jsonencode(s), 'char');
end
