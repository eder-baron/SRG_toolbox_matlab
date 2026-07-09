function srg_plot_gauss(TransferFcn, gplus, gminus, color, estpoints, options)
%SRG_PLOT_GAUSS  Plot the SRG in the Gauss (complex) plane.
%
%   SRG_PLOT_GAUSS(TransferFcn, gplus, gminus, color, estpoints)
%
%   Name-Value Arguments:
%       Theme     - Style theme (default: "default")
%       FaceAlpha - Fill transparency (default: 0.5)
%       Fill      - Fill the interior of each slice (default: false)
%                   Useful for MIMO systems to see the SRG as a filled region.
%
%   SISO detection:
%       Checks whether every cell in gplus contains a single unique complex
%       value, regardless of TransferFcn type. This handles all cases:
%       scalar tf, ss, and even numeric scalars passed via srg_compute.
%
%   See also SRG_COMPUTE, SRG_PLOT_BELTRAMI, SRG_PLOT_COMPARE

    arguments
        TransferFcn
        gplus       cell
        gminus      cell
        color
        estpoints   (1,1) double
        options.Theme     string  = "default"
        options.FaceAlpha (1,1) double = 0.5
        options.Fill      (1,1) logical = false
    end

    s = srg_style('Theme', options.Theme);
    c = resolve_color(color, s);
    c_pale = c + (1 - c) * 0.45;

    hold on

    n = min(estpoints, numel(gplus));

    % =====================================================================
    % Route based on data shape — data-driven, not TransferFcn-type-driven
    % =====================================================================
    if isa(TransferFcn, 'double')
        plot_constant(gplus, gminus, c, s, options);

    elseif is_siso_cell(gplus, n)
        plot_siso(gplus, gminus, n, c, c_pale, s, options);

    else
        plot_mimo(gplus, gminus, n, c, c_pale, s, options);
    end

    srg_apply_style(s);

end


% =========================================================================
%  PLOTTING BRANCHES
% =========================================================================

function plot_constant(gplus, gminus, c, s, options)
%PLOT_CONSTANT  Static matrix: single boundary or marker.
    if is_single_point(gplus{1})
        % Scalar constant — just markers
        plot(real(gplus{1}(1)),  imag(gplus{1}(1)),  'x', ...
             'Color', c, 'LineWidth', s.linewidth, 'MarkerSize', 10)
        plot(real(gminus{1}(1)), imag(gminus{1}(1)), 'x', ...
             'Color', c, 'LineWidth', s.linewidth, 'MarkerSize', 10)
    else
        % Matrix constant — draw (and optionally fill) the single boundary
        if options.Fill
            fill_slice(gplus{1}, c, options.FaceAlpha);
        end
        plot(real(gplus{1}),  imag(gplus{1}),  'Color', c, 'LineWidth', s.linewidth)
        plot(real(gminus{1}), imag(gminus{1}), 'Color', c, 'LineWidth', s.linewidth)
    end
end


function plot_siso(gplus, gminus, n, c, c_pale, s, options)
%PLOT_SISO  SISO: collect one point per frequency and draw trajectories.
%
%   The SRG of a SISO system is the region swept between the gplus
%   (upper) and gminus (lower = conj(gplus)) traces as frequency varies.

    gp_trace = zeros(n, 1);
    gm_trace = zeros(n, 1);
    for ii = 1:n
        gp_trace(ii) = first_finite(gplus{ii});
        gm_trace(ii) = first_finite(gminus{ii});
    end

    % Remove any frequencies that produced NaN
    valid     = isfinite(gp_trace) & isfinite(gm_trace);
    gp_trace  = gp_trace(valid);
    gm_trace  = gm_trace(valid);

    if isempty(gp_trace), return; end

    if options.Fill
        % Closed loop: upper sweep forward, lower sweep backward
        boundary = [gp_trace; flipud(gm_trace)];
        ps = polyshape(real(boundary), imag(boundary), 'Simplify', false);
        if ps.NumRegions > 0
            plot(ps, 'FaceColor', c_pale, 'FaceAlpha', options.FaceAlpha, ...
                 'EdgeColor', c, 'EdgeAlpha', 0.7, 'LineWidth', s.linewidth);
        end
    end

    % Upper (gplus) and lower (gminus) boundary curves — two lines only
    plot(real(gp_trace), imag(gp_trace), '-', 'Color', c, 'LineWidth', s.linewidth)
    plot(real(gm_trace), imag(gm_trace), '-', 'Color', c, 'LineWidth', s.linewidth)
end


function plot_mimo(gplus, gminus, n, c, c_pale, s, options)
%PLOT_MIMO  MIMO: per-frequency curves, optionally with filled slices.
    for ii = 1:n
        if isempty(gplus{ii}), continue; end

        if options.Fill
            fill_slice(gplus{ii}, c_pale, options.FaceAlpha);
        end

        plot(real(gplus{ii}), imag(gplus{ii}), ...
             'LineWidth', s.linewidth, 'Color', c)

        if ii <= numel(gminus) && ~isempty(gminus{ii})
            plot(real(gminus{ii}), imag(gminus{ii}), ...
                 'LineWidth', s.linewidth, 'Color', c)
        end
    end
end


% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function tf = is_siso_cell(gp, n)
%IS_SISO_CELL  True when every non-empty cell contains one unique point.
%   Works regardless of TransferFcn class. Handles the srg_compute case
%   where field_of_values() returns [e; e] for scalar inputs (length 2, all identical).
    tf  = true;
    tol = 1e-8;
    for ii = 1:n
        c = gp{ii}(isfinite(gp{ii}));
        if isempty(c), continue; end
        if max(abs(c - c(1))) > tol * max(abs(c(1)), 1)
            tf = false;
            return;
        end
    end
end


function tf = is_single_point(curve)
%IS_SINGLE_POINT  True when all finite values in curve are the same point.
    curve = curve(isfinite(curve));
    if isempty(curve), tf = true; return; end
    tol = max(1e-10, max(abs(curve)) * 1e-8);
    tf  = max(abs(curve - curve(1))) < tol;
end


function v = first_finite(curve)
%FIRST_FINITE  Return the first finite value in curve, or NaN.
    curve = curve(isfinite(curve));
    if isempty(curve), v = NaN; else, v = curve(1); end
end


function fill_slice(gp, c_pale, face_alpha)
%FILL_SLICE  Fill the region bounded by gplus (and its conjugate mirror).
%   For the Gauss plane, the SRG at one frequency is symmetric about the
%   real axis: the boundary is gp ∪ conj(flip(gp)).
    boundary = gp(:);
    boundary = boundary(isfinite(boundary));
    if length(boundary) < 3, return; end

    tol  = max(1e-15, max(abs(boundary)) * 1e-10);
    keep = [true; abs(diff(boundary)) > tol];
    boundary = boundary(keep);
    if length(boundary) < 3, return; end

    % Close with conjugate to get the full symmetric region
    boundary_full = [boundary; conj(flipud(boundary))];

    ps = polyshape(real(boundary_full), imag(boundary_full), 'Simplify', false);
    if ps.NumRegions > 0
        plot(ps, 'FaceColor', c_pale, 'FaceAlpha', face_alpha, ...
             'EdgeColor', 'none');
    end
end


function c = resolve_color(color, s)
    if isnumeric(color) && isscalar(color)
        c = s.color(color);
    elseif isnumeric(color) && numel(color) == 3
        c = color(:)';
    elseif isstring(color) || ischar(color)
        c = color;
    else
        c = s.color(1);
    end
end