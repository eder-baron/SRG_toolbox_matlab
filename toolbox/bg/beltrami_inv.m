function [gplus, gminus] = beltrami_inv(T)
%BELTRAMI_INV  Inverse Beltrami-Klein mapping to the Gauss plane.
%
%   [GPLUS, GMINUS] = BELTRAMI_INV(T) maps points from the Beltrami-Klein
%   disk back to the Gauss (complex) plane using:
%       g+(z) = (Im(z) + i*sqrt(1 - |z|^2)) / (Re(z) - 1)
%       g-(z) = (Im(z) - i*sqrt(1 - |z|^2)) / (Re(z) - 1)
%
%   Input:
%       T   - Matrix of complex values in the Beltrami-Klein disk
%
%   Outputs:
%       gplus  - Upper branch of the inverse mapping
%       gminus - Lower branch of the inverse mapping
%
%   See also BELTRAMI_MAP_SCALAR, BELTRAMI_MAP_MATRIX

gplus  = zeros(size(T));
gminus = zeros(size(T));
[a, b] = size(T);
for ii = 1:a
    for jj = 1:b
        gplus(ii,jj)  = (imag(T(ii,jj)) + 1i*sqrt(1 - norm(T(ii,jj))^2)) ...
                        / (real(T(ii,jj)) - 1);
        gminus(ii,jj) = (imag(T(ii,jj)) - 1i*sqrt(1 - norm(T(ii,jj))^2)) ...
                        / (real(T(ii,jj)) - 1);
    end
end

end
