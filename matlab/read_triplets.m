function M = read_triplets(path)
%READ_TRIPLETS Read a sparse matrix written by phcpf_mesher: header "rows cols nnz", then 0-based
%   "i j value" lines.
fid = fopen(path, 'r');
if fid < 0, error('read_triplets: cannot open %s', path); end
hdr = fscanf(fid, '%d', 3);
T = fscanf(fid, '%f', [3, hdr(3)])';
fclose(fid);
if isempty(T), T = zeros(0, 3); end
M = sparse(T(:, 1) + 1, T(:, 2) + 1, T(:, 3), hdr(1), hdr(2));
end
