function [srg, srg_upper] = srg_hard(G, alphas, phis, frequencies, options)
%SRG_HARD  Compute the hard SRG boundary of an LTI system.
%
%   [SRG, SRG_UPPER] = SRG_HARD(G, ALPHAS, PHIS, FREQUENCIES) computes
%   the boundary of the hard Scaled Relative Graph.
%
%   The hard SRG differs from the soft SRG in that it accounts for:
%       - Unstable systems: max gain -> infval (large finite value)
%       - Non-minimum-phase zeros: min gain -> 0
%
%   Inputs:
%       G           - Transfer function (tf/ss object)
%       alphas      - Vector of real-valued base points
%       phis        - Vector of angles in [0, pi] for boundary resolution
%       frequencies - Vector of frequencies in rad/s
%
%   Name-Value Arguments:
%       InfVal - Finite replacement for infinite gain (default: 10)
%
%   Outputs:
%       srg       - Complex vector: full SRG boundary
%       srg_upper - Complex vector: upper half-plane boundary only
%
%   Example:
%       G = tf([1 -1], [1 2 1]);  % non-minimum phase
%       alphas = linspace(-2, 2, 50);
%       phis = linspace(0, pi, 200);
%       freqs = logspace(-3, 3, 500);
%       [srg_boundary, ~] = srg_hard(G, alphas, phis, freqs);
%
%       figure; fill(real(srg_boundary), imag(srg_boundary), 'k', ...
%                    'FaceAlpha', 0.4, 'EdgeColor', 'none');
%       axis equal; grid on; title('Hard SRG');
%
%   Based on SrgTools.jl by Julius P. J. Krebbekx
%   (BSD 3-Clause License).

    arguments
        G
        alphas      (1,:) double
        phis        (1,:) double
        frequencies (1,:) double
        options.InfVal (1,1) double {mustBePositive} = 10
    end

    % Step 1: Compute min/max gain circles (hard version)
    [min_radii, max_radii] = srg_circles(G, alphas, frequencies, options.InfVal);

    % Step 2: Find a real point z0 inside all max-gain circles
    [a, b] = circle_interval(alphas, max_radii);
    z0 = (a + b) / 2;

    % Step 3: Intersect max-gain circles
    srg_max_radius = max_boundary(z0, alphas, max_radii, phis);

    % Step 4: Carve out min-gain circles
    srg_upper = carve_min(z0, srg_max_radius, phis, alphas, min_radii);

    % Step 5: Mirror to get full SRG
    srg = [srg_upper; conj(flipud(srg_upper))];

end

% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function [min_radii, max_radii] = srg_circles(G, alphas, frequencies, infval)
%COMPUTE_SRG_CIRCLES_HARD  Hard SRG circles with stability and zero checks.

    [q, p] = size(G);

    if q < p
        G_padded = [G; tf(zeros(p - q, p))];
    else
        G_padded = G;
    end

    n_rows = max(q, p);
    Id = tf(eye(n_rows, p));

    N = length(alphas);
    min_radii = zeros(N, 1);
    max_radii = zeros(N, 1);

    for ii = 1:N
        G_a = G_padded - alphas(ii) * Id;
        [min_radii(ii), max_radii(ii)] = gain_range(G_a, frequencies, infval);
    end

end

% -------------------------------------------------------------------------

function [sigma_min, sigma_max] = gain_range(G, frequencies, infval)
%COMPUTE_MIN_MAX_GAIN_HARD  Min/max gains with stability and zero checks.

    G_jw = freqresp(G, 1i*frequencies(1));
    sigmas = svd(G_jw);
    sigma_min = min(sigmas);
    sigma_max = max(sigmas);

    for ii = 2:length(frequencies)
        G_jw = freqresp(G, 1i*frequencies(ii));
        sigmas = svd(G_jw);
        sigma_max = max(sigma_max, max(sigmas));
        sigma_min = min(sigma_min, min(sigmas));
    end

    if ~isstable(G)
        sigma_max = infval;
    end

    z = tzero(G);
    if ~isempty(z)
        if max(real(z)) >= 0
            sigma_min = 0;
        end
    end

end

% -------------------------------------------------------------------------

function [a, b] = circle_interval(alphas, max_radii)
%FIND_CIRCLE_INTERVAL  Real interval contained in all max-gain circles.

    a = alphas(1) - max_radii(1);
    b = alphas(1) + max_radii(1);

    for ii = 2:length(alphas)
        a = max(a, alphas(ii) - max_radii(ii));
        b = min(b, alphas(ii) + max_radii(ii));
    end

end

% -------------------------------------------------------------------------

function srg_max_radius = max_boundary(z0, alphas, max_radii, phis)
%BOUNDARY_MAX_GAIN_INTERSECTIONS  Outer radius from max-gain circle intersections.

    n = length(phis);
    m = length(alphas);
    srg_max_radius = zeros(n, 1);
    radii = zeros(m, 1);

    for ii = 1:n
        for jj = 1:m
            varphi = varphi_center(z0, alphas(jj), max_radii(jj), phis(ii));
            radii(jj) = abs(alphas(jj) + max_radii(jj)*exp(1i*varphi) - z0);
        end
        srg_max_radius(ii) = min(radii);
    end

end

% -------------------------------------------------------------------------

function varphi = varphi_center(z0, alpha, r_alpha, phi)
%COMPUTE_VARPHI_SRG_CENTER  Angle conversion between polar decompositions.

    assert(phi >= 0 && phi <= pi, 'phi must be in [0, pi]');

    if phi == pi/2
        varphi = acos(-(alpha - z0) / r_alpha);
        return;
    end

    if phi > pi/2
        phi = phi - pi;
    end

    tanphi2 = tan(phi)^2;
    a = r_alpha^2 * tanphi2 + r_alpha^2;
    b = 2 * r_alpha * (alpha - z0) * tanphi2;
    c = (alpha - z0)^2 * tanphi2 - r_alpha^2;

    disc = b^2 - 4*a*c;
    x1 = (-b - sqrt(disc)) / (2*a);
    x2 = (-b + sqrt(disc)) / (2*a);

    if abs(x1) > 1 && abs(x2) <= 1
        varphi = acos(x2);
    elseif abs(x2) > 1 && abs(x1) <= 1
        varphi = acos(x1);
    else
        err1 = abs(tan(phi)*(alpha - z0 + r_alpha*x1) - r_alpha*sqrt(1 - x1^2));
        err2 = abs(tan(phi)*(alpha - z0 + r_alpha*x2) - r_alpha*sqrt(1 - x2^2));
        if err1 < err2
            varphi = acos(x1);
        else
            varphi = acos(x2);
        end
    end

end

% -------------------------------------------------------------------------

function srg_outer = carve_min(z0, srg_max_radius, phis, alphas, min_radii)
%REMOVE_MIN_CIRCLES  Carve min-gain circles from the max-gain boundary.

    srg_outer = z0 + srg_max_radius(:) .* exp(1i * phis(:));

    a = min(real(srg_outer));
    b = max(real(srg_outer));

    real_line = linspace(a, b, length(phis))';
    srg_outer = [srg_outer; real_line];

    % First pass: centers outside [a, b]
    for ii = 1:length(alphas)
        if alphas(ii) < a || alphas(ii) > b
            for jj = 1:length(srg_outer)
                z = srg_outer(jj) - alphas(ii);
                varphi = atan2(imag(z), real(z));
                if abs(z) < min_radii(ii)
                    srg_outer(jj) = alphas(ii) + min_radii(ii) * exp(1i*varphi);
                end
            end
        end
    end

    % Second pass: centers inside [a, b]
    for ii = 1:length(alphas)
        if alphas(ii) >= a && alphas(ii) <= b
            for jj = 1:length(srg_outer)
                z = srg_outer(jj) - alphas(ii);
                if abs(z) < min_radii(ii)
                    re_z = real(z);
                    srg_outer(jj) = alphas(ii) + re_z + ...
                        1i * sqrt(1 - (re_z/min_radii(ii))^2) * min_radii(ii);
                end
            end
        end
    end

end