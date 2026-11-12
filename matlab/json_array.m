function s = json_array(x, dtype)
%JSON_ARRAY Pack a numeric matrix for the web viewer as {dtype, shape, b64}.
%   Rows are flattened in row-major order (so an n x 3 vertex matrix becomes x1 y1 z1 x2 ...).
%   dtype: 'f32' (default), 'i32' (0-based indices are the caller's job) or 'u8'.
if nargin < 2, dtype = 'f32'; end
if isempty(x), x = zeros(0, 1); end
if isvector(x) && ~ismatrix(x), x = x(:); end
shape = size(x);
flat = reshape(x.', 1, []);
switch dtype
    case 'f32', bytes = typecast(single(flat), 'uint8');
    case 'i32', bytes = typecast(int32(flat), 'uint8');
    case 'u8',  bytes = uint8(flat);
    otherwise, error('json_array: unknown dtype %s', dtype);
end
s = struct('dtype', dtype, 'shape', shape, 'b64', matlab.net.base64encode(bytes));
end
