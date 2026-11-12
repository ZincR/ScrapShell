function [V, F] = read_poly_off(path)
%READ_POLY_OFF Read a polygonal OFF file. F is a cell array of 1-based vertex index rows.
txt = fileread(path);
lines = regexp(txt, '\r?\n', 'split');
lines = lines(~cellfun(@isempty, strtrim(lines)));
lines = lines(~startsWith(strtrim(lines), '#'));
k = 1;
if startsWith(strtrim(lines{k}), 'OFF')
    rest = strtrim(erase(lines{k}, 'OFF'));
    if isempty(rest), k = k + 1; counts = sscanf(lines{k}, '%d'); else, counts = sscanf(rest, '%d'); end
else
    error('read_poly_off: %s is not an OFF file', path);
end
nv = counts(1); nf = counts(2);
V = zeros(nv, 3);
for i = 1:nv
    V(i, :) = sscanf(lines{k + i}, '%f', 3)';
end
F = cell(nf, 1);
for i = 1:nf
    row = sscanf(lines{k + nv + i}, '%d')';
    F{i} = row(2:row(1) + 1) + 1;
end
end
