function Wrefined = srg_refine_chords(W, points, tol)
%SRG_REFINE_CHORDS  Subdivide long chords in a Beltrami-Klein boundary.
%
%   WREFINED = SRG_REFINE_CHORDS(W, POINTS, TOL) inserts additional
%   points along any segment of the boundary W whose real- or
%   imaginary-part step exceeds TOL, by linearly interpolating POINTS
%   points along the chord connecting the two endpoints.
%
%   This is required before calling BELTRAMI_INV. FIELD_OF_VALUES
%   returns the *exact* Beltrami-Klein boundary for normal matrices
%   (e.g. real-diagonal operators), since for a normal matrix the field
%   of values is precisely the convex hull of its eigenvalues, and
%   straight polygon edges are correct there. BELTRAMI_INV, however, is
%   a strongly nonlinear map
%       g(z) = (Im(z) +- i*sqrt(1-|z|^2)) / (Re(z) - 1)
%   so a straight BK-plane chord does not map to a straight Gauss-plane
%   chord. Passing an under-sampled boundary (e.g. only the handful of
%   vertex points FIELD_OF_VALUES returns for a normal matrix) through
%   BELTRAMI_INV and connecting the *mapped* outputs with straight
%   lines produces a badly distorted, collapsed-looking curve, even
%   though the BK-plane data was exact. Densifying long chords here,
%   entirely within the BK disk (a linear space), fixes this without
%   touching FIELD_OF_VALUES or BELTRAMI_INV.
%
%   Inputs:
%       W      - Column vector of complex boundary points (BK-plane)
%       points - Number of interpolated points to insert along a chord
%                that exceeds TOL (same convention as the angular
%                resolution passed to FIELD_OF_VALUES)
%       tol    - Maximum allowed real- or imaginary-part step between
%                consecutive boundary points before subdivision
%
%   Output:
%       Wrefined - Boundary with long chords subdivided. If no chord
%                  exceeds TOL, Wrefined is identical to W.
%
%   See also FIELD_OF_VALUES, BELTRAMI_INV

    W = W(:);

    Wr = refine_by_component(W, points, tol, 'real');
    Wrefined = refine_by_component(Wr, points, tol, 'imag');

end

%--------------------------------------------------------------------------
function Wout = refine_by_component(W, points, tol, component)
%REFINE_BY_COMPONENT  Subdivide chords whose step in REAL or IMAG part
%   exceeds TOL, walking backwards so inserted indices don't disturb the
%   remaining (not-yet-processed) segments.

    Wout = W;
    switch component
        case 'real'
            comp = real(W);
        case 'imag'
            comp = imag(W);
    end

    for ii = length(W):-1:2
        if abs(comp(ii) - comp(ii-1)) > tol
            x1 = real(W(ii-1)); x2 = real(W(ii));
            y1 = imag(W(ii-1)); y2 = imag(W(ii));

            % Two points determine the chord exactly -- no fit needed.
            % Parameterize along whichever axis has the larger spread,
            % to avoid dividing by a near-zero run (a near-vertical
            % chord parameterized by x, or a near-horizontal one by y).
            if abs(x2 - x1) >= abs(y2 - y1)
                xi = linspace(x1, x2, points)';
                yi = y1 + (xi - x1) * (y2 - y1) / (x2 - x1);
            else
                yi = linspace(y1, y2, points)';
                xi = x1 + (yi - y1) * (x2 - x1) / (y2 - y1);
            end

            lines = xi + 1i*yi;
            Wout = [Wout(1:ii-1); lines; Wout(ii:end)];
        end
    end

end
