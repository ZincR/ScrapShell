function write_poly_off(path, V, F)
%WRITE_POLY_OFF Write a polygonal OFF file. F is a cell array of 1-based index rows.
fid = fopen(path, 'w');
if fid < 0, error('write_poly_off: cannot open %s', path); end
cleaner = onCleanup(@() fclose(fid));
fprintf(fid, 'OFF\n%d %d 0\n', size(V, 1), numel(F));
fprintf(fid, '%.17g %.17g %.17g\n', V');
for i = 1:numel(F)
    f = F{i};
    fprintf(fid, '%d', numel(f));
    fprintf(fid, ' %d', f - 1);
    fprintf(fid, '\n');
end
end
