function map = beltrami_map_scalar(T)
%BELTRAMI_MAP_SCALAR  Beltrami-Klein mapping for scalar eigenvalues.
%
%   MAP = BELTRAMI_MAP_SCALAR(T) maps complex values T into the
%   Beltrami-Klein disk using the scalar formula:
%       map(i,j) = (conj(T)-i)*(T-i) / (1 + conj(T)*T)
%
%   Input:
%       T   - Matrix of complex values (typically eigenvalues)
%
%   Output:
%       map - Mapped values in the Beltrami-Klein disk
%
%   See also BELTRAMI_MAP_MATRIX, BELTRAMI_INV

map = zeros(size(T));
[a, b] = size(T);
for ii = 1:a
    for jj = 1:b
        map(ii,jj) = (conj(T(ii,jj)) - 1i) * (T(ii,jj) - 1i) ...
                     / (1 + conj(T(ii,jj)) * T(ii,jj));
    end
end

end
