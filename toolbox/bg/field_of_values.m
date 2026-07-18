function [f, e] = field_of_values(B, nk, thmax, noplot)
%FIELD_OF_VALUES  Field of values (numerical range) of a matrix.
%
%   FIELD_OF_VALUES(A, NK, THMAX) evaluates and plots the field of values
%   of the NK largest leading principal submatrices of A, using THMAX
%   equally spaced angles in the complex plane.
%   The defaults are NK = 1 and THMAX = 16.
%   (For publication quality, set THMAX higher, e.g. 32.)
%   The eigenvalues of A are displayed as 'x'.
%
%   [F, E] = FIELD_OF_VALUES(A, NK, THMAX, 1) suppresses the plot and
%   returns the field of values boundary in F, with A's eigenvalues in E.
%   NORM(F,INF) approximates the numerical radius.
%
%   THMAX = 1 (e.g. from srg_compute's SISO convention of points = 1)
%   requests a single angle sample. The rotation loop below is
%   parameterized by thmax-1 = 0 angle steps, which would require
%   dividing by zero; in that case F reduces to A's eigenvalues
%   directly (E), which is checked immediately after E is computed,
%   before the rotation loop runs.
%
%   Based on an original routine by A. Ruhe.
%
%   References:
%       R. A. Horn and C. R. Johnson, Topics in Matrix Analysis, Cambridge
%            University Press, 1991; sec. 1.5.
%       A. S. Householder, The Theory of Matrices in Numerical Analysis,
%            Blaisdell, New York, 1964; sec. 3.3.
%       C. R. Johnson, Numerical determination of the field of values of a
%            general complex matrix, SIAM J. Numer. Anal., 15 (1978),
%            pp. 595-602.

if nargin < 2 || isempty(nk),    nk = 1;    end
if nargin < 3 || isempty(thmax), thmax = 16; end
thmax = thmax - 1;  % Code below uses thmax + 1 angles.

iu = 1i;
[n, p] = size(B);
if n ~= p, error('field_of_values:notSquare', 'Matrix must be square.'); end
f = [];
e = eig(B);

if thmax == 0
    % A single angle sample was requested (THMAX input = 1). The
    % rotation loop below computes th = ii/thmax*pi, which is 0/0 here
    % -- so skip it entirely and use the eigenvalues directly, matching
    % the fallback that used to sit (unreachably, for this case) after
    % the loop.
    f = e;

elseif isequal(B, B')
    % Hermitian: field of values is the real interval [min(e), max(e)].
    f = [min(e) max(e)];

elseif isequal(B, -B')
    % Skew-Hermitian: field of values is a purely imaginary interval.
    e = imag(e);
    f = [min(e) max(e)];
    e = iu*e; f = iu*f;

else
    z = zeros(2*thmax+1, 1);
    for m = 1:nk
        ns = n + 1 - m;
        A = B(1:ns, 1:ns);

        for ii = 0:thmax
            th = ii / thmax * pi;
            Ath = exp(iu*th) * A;             % Rotate A through angle th.
            H = 0.5 * (Ath + Ath');           % Hermitian part of rotated A.
            [X, D] = eig(H);
            [~, k] = sort(real(diag(D)));
            z(1+ii)       = rq(A, X(:,k(1)));   % Smallest real part
            z(1+ii+thmax) = rq(A, X(:,k(ns)));  % Largest real part
        end

        f = [f; z]; %#ok<AGROW>
    end
    % Ensure boundary is closed (needed for orthogonal matrices).
    f = [f; f(1,:)];
end

% Plot if no output suppression argument given
if nargin < 4
    % Compute balanced axis limits from the data (replaces cpltaxes)
    ax = compute_axis_limits(f);

    plot(real(f), imag(f))
    axis(ax);
    axis('square');

    hold on
    plot(real(e), imag(e), 'x')
    hold off
end

end

%--------------------------------------------------------------------------
function z = rq(A, x)
%RQ  Rayleigh quotient: x'*A*x / (x'*x).
    z = x'*A*x / (x'*x);
end

%--------------------------------------------------------------------------
function ax = compute_axis_limits(f)
%COMPUTE_AXIS_LIMITS  Balanced square axis limits for complex data.
%   Replaces the external cpltaxes dependency.
    rp = real(f(:));
    ip = imag(f(:));
    margin = 0.1;

    xmin = min(rp); xmax = max(rp);
    ymin = min(ip); ymax = max(ip);

    dx = xmax - xmin;
    dy = ymax - ymin;

    if dx == 0 && dy == 0
        ax = [xmin-1 xmax+1 ymin-1 ymax+1];
        return;
    end

    % Make square and add margin
    span = max(dx, dy);
    cx = (xmin + xmax) / 2;
    cy = (ymin + ymax) / 2;
    half = span / 2 * (1 + margin);

    ax = [cx-half cx+half cy-half cy+half];
end

