function map = beltrami_map_matrix(T)
%BELTRAMI_MAP_MATRIX  Beltrami-Klein mapping for a matrix operator.
%
%   MAP = BELTRAMI_MAP_MATRIX(T) maps the matrix operator T into the
%   Beltrami-Klein disk using the matrix formula:
%       map = (I + T'*T)^(-1/2) * (T' - iI) * (T - iI) * (I + T'*T)^(-1/2)
%
%   Input:
%       T   - Square complex matrix
%
%   Output:
%       map - Mapped matrix in the Beltrami-Klein disk
%
%   See also BELTRAMI_MAP_SCALAR, BELTRAMI_INV

I = eye(length(T));
map = (I + T'*T)^(-1/2) * (T' - 1i*I) * (T - 1i*I) * (I + T'*T)^(-1/2);

end
