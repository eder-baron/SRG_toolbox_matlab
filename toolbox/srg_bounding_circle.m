function [center, radius] = srg_bounding_circle(srg, options)
%SRG_BOUNDING_CIRCLE  Smallest circle enclosing all SRG points.
%
%   [CENTER, RADIUS] = SRG_BOUNDING_CIRCLE(SRG) computes the minimum
%   enclosing circle (smallest circle containing all points) for a set
%   of complex-valued SRG boundary points.
%
%   Since SRGs are symmetric about the real axis, the center is
%   constrained to the real line.
%
%   Name-Value Arguments:
%       Plot - If true, plot the circle over the current axes (default: false)
%
%   Inputs:
%       srg - Complex vector of SRG boundary points
%
%   Outputs:
%       center - Complex scalar (real-valued) center of the circle
%       radius - Radius of the minimum enclosing circle
%
%   Example:
%       G = tf([1], [1 1]);
%       [~, ~, gplus, ~, ~] = srg_compute(G, -2, 3, 200, 32);
%       [sg_plus, ~, ~] = srg_scaled_graph([], 'InputDomain', 'Gauss', ...
%                             'GPlus', gplus, 'GMinus', gplus);
%       [c, r] = srg_bounding_circle(sg_plus{1}(:), 'Plot', true);
%       fprintf('Center = %.4f, Radius = %.4f\n', c, r);
%
%   See also SRG_HARD, SRG_QSR_MULTIPLIER

    arguments
        srg (:,1) double
        options.Plot (1,1) logical = false
    end

    % Work with upper half-plane only (symmetry guarantees same result)
    % but keep all points for the actual enclosing computation
    x = real(srg);
    y = imag(srg);

    % Since SRG is symmetric about the real axis, the optimal center
    % lies on the real axis. We solve:
    %   min_{c real} max_i |srg(i) - c|
    %
    % This is a 1D minimax problem. For a point z = x + iy and real center c:
    %   |z - c|^2 = (x - c)^2 + y^2
    %
    % The optimal c minimizes the maximum of these distances.
    % Use fminbnd on the real extent.

    c_lo = min(x) - max(abs(y));
    c_hi = max(x) + max(abs(y));

    % Objective: maximum distance from real center c to any point
    max_dist = @(c) max(sqrt((x - c).^2 + y.^2));

    % Golden section search (fminbnd)
    [c_opt, r_opt] = fminbnd(max_dist, c_lo, c_hi, ...
                              optimset('TolX', 1e-12, 'Display', 'off'));

    center = c_opt;
    radius = r_opt;

    % Plot if requested
    if options.Plot
        theta = linspace(0, 2*pi, 500);
        z_circle = center + radius * exp(1i * theta);

        hold on
        plot(real(z_circle), imag(z_circle), 'r-', 'LineWidth', 1.5);
        plot(center, 0, 'r+', 'MarkerSize', 10, 'LineWidth', 1.5);
        hold off
    end

end