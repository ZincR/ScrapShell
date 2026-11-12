function S = fast_sparse(varargin)
%FAST_SPARSE Drop-in for the Windows-only fast_sparse.mexw64 shipped with PH-CPF.
%   Same signature as SPARSE (duplicates are summed), so we simply forward.
S = sparse(varargin{:});
end
