function [f, e] = field_of_values(B, nk, thmax, noplot)
%FIELD_OF_VALUES  Field of values (numerical range) of a matrix.
%
%   [F, E] = FIELD_OF_VALUES(B, NK, THMAX, 1) computes the field of values
%   boundary of the NK largest leading principal submatrices of B using
%   THMAX equally spaced angles, and returns the boundary in F and the
%   eigenvalues of B in E.  Always pass the fourth argument as 1 to
%   suppress plotting; the plotting branch requires cpltaxes from the
%   Matrix Computation Toolbox and has been removed from this copy.
%
%   Defaults: NK = 1, THMAX = 16.
%   For smoother SRG boundaries use THMAX = 64 or higher.
%
%   NORM(F, INF) approximates the numerical radius
%       max { |z| : z in field_of_values(B) }.
%
%   Theory:
%       field_of_values(A) is the set of all Rayleigh quotients x'Ax/(x'x).
%       It is a convex set containing the eigenvalues of A.  When A is
%       normal, field_of_values(A) equals the convex hull of the eigenvalues.
%           z = x'Ax/(x'x)  =>  Re(z) = x'Hx/(x'x),  H = (A+A')/2
%       so  min(eig(H)) <= Re(z) <= max(eig(H)), with equality for the
%       corresponding eigenvectors of H.
%
%   Based on an original routine by A. Ruhe (via N. J. Higham,
%   Matrix Computation Toolbox).
%
%   References:
%       R. A. Horn and C. R. Johnson, Topics in Matrix Analysis,
%           Cambridge University Press, 1991, sec. 1.5.
%       C. R. Johnson, Numerical determination of the field of values of a
%           general complex matrix, SIAM J. Numer. Anal., 15 (1978), 595-602.
%
%   See also BELTRAMI_MAP_MATRIX, SRG_COMPUTE

if nargin < 2 || isempty(nk),    nk    = 1;  end
if nargin < 3 || isempty(thmax), thmax = 16; end
if nargin < 4,                   noplot = false; end

if ~noplot
    error('field_of_values:noPlot', ...
        ['The plotting branch of field_of_values requires cpltaxes from\n' ...
         'the Matrix Computation Toolbox, which is not bundled with the\n' ...
         'SRG Toolbox.  Pass a fourth argument (e.g. 1) to suppress\n' ...
         'plotting and return [f, e] instead.']);
end

thmax = thmax - 1;  % code below uses thmax+1 angles

iu = 1i;
[n, p] = size(B);
if n ~= p, error('field_of_values:nonSquare', 'Matrix must be square.'); end

f = [];
z = zeros(2*thmax + 1, 1);
e = eig(B);

% Shortcuts for Hermitian and skew-Hermitian matrices
if isequal(B, B')
    f = [min(e); max(e)];

elseif isequal(B, -B')
    e = imag(e);
    f = [min(e); max(e)];
    e = iu * e;
    f = iu * f;

else
    for m = 1:nk
        ns  = n + 1 - m;
        A   = B(1:ns, 1:ns);

        for i = 0:thmax
            th  = i / thmax * pi;
            Ath = exp(iu * th) * A;         % rotate A through angle th
            H   = 0.5 * (Ath + Ath');       % Hermitian part
            [X, D] = eig(H);
            [~, k]  = sort(real(diag(D)));
            z(1 + i)        = rayleigh_quotient(A, X(:, k(1)));   % smallest Re
            z(1 + i + thmax) = rayleigh_quotient(A, X(:, k(ns))); % largest  Re
        end

        f = [f; z]; %#ok<AGROW>
    end

    % Close the boundary (needed for orthogonal matrices)
    f = [f; f(1, :)];
end

if thmax == 0
    f = e;
end

end % field_of_values


% -------------------------------------------------------------------------
function z = rayleigh_quotient(A, x)
%RAYLEIGH_QUOTIENT  z = x'*A*x / (x'*x).
z = (x' * A * x) / (x' * x);
end
