function [mtch, S] = write_cpf_matching(cpf, prefix)
%WRITE_CPF_MATCHING Export the CPF matching and singularities for tools/phcpf_mesher, as the official
%   pipeline passed them to MeshCutter_bin (CPF.uv2mesh -> compute_matching_and_sing -> cut_and_comb).
%   <prefix>_matching.txt: "v0 v1 fi fj m" per inner edge (0-based), vector k of face fi matches
%   vector k+m of face fj. <prefix>_singin.txt: "v index" per singular vertex (0-based).
M = cpf.M;
[mtch, S] = cpf.compute_matching_and_sing(cpf.U, cpf.V);
e = M.inner_edges(:);
rows = [M.edges(e, 1) - 1, M.edges(e, 2) - 1, abs(M.e2t(e, 1)) - 1, abs(M.e2t(e, 2)) - 1, mtch(e)];
writematrix(rows, [prefix '_matching.txt'], 'Delimiter', ' ', 'FileType', 'text');
sv = find(S);
writematrix([sv - 1, S(sv)], [prefix '_singin.txt'], 'Delimiter', ' ', 'FileType', 'text');
end
