function srg_plot_3d(rawdata, gplus1, fmin, fmax, color, options)
%SRG_PLOT_3D  3D wire/filled plot of SRG boundaries across frequencies.
%
%   SRG_PLOT_3D(rawdata, gplus1, fmin, fmax, color)
%
%   For SISO systems (scalar gplus at each frequency), collects all
%   values into a single 3D trace. For MIMO systems, plots each
%   frequency slice as a separate curve or filled patch.
%   For constant matrices (single frequency), plots with visible markers.
%
%   Inputs:
%       rawdata - Table with .freq field
%       gplus1  - Cell array of SRG boundary data per frequency
%       fmin    - Minimum frequency to display (0 = first available)
%       fmax    - Maximum frequency to display (0 = last available)
%       color   - Color index (1-7) or RGB triplet
%
%   Name-Value Arguments:
%       Theme         - Style theme (default: "default")
%       Mirror        - Plot conjugate (lower) half (default: true)
%       gminus        - Cell array of lower SRG boundaries for SISO (default: {})
%       Fill          - Fill each frequency slice as a solid patch (default: false)
%                       Applies to MIMO only; silently ignored (with warning) for SISO.
%       FaceAlpha     - Transparency of filled patches, in [0,1] (default: 0.3)
%       ShowBoundary  - Overlay boundary wire on top of fill (default: true)
%                       Only meaningful when Fill = true; wire is always drawn
%                       when Fill = false.
%       BoundaryColor - Color of the boundary wire when Fill = true.
%                       Accepts an RGB triplet or a color index (1-7).
%                       Default: same color as the fill (i.e., color argument).
%       FillTightness - Controls how tightly the filled polygon follows the
%                       boundary data, in [0, 1] (default: 0).
%                         0 → convex hull of the slice points (always convex).
%                         1 → exact boundary shape via polyshape on the
%                             ordered curve (can be non-convex).
%                       Values between 0 and 1 are rounded to the nearest
%                       extreme so the semantics stay well-defined.
%
%   See also SRG_COMPUTE, SRG_PLOT_3D_SURFACE, SRG_PLOT_3D_BELTRAMI

    arguments
        rawdata
        gplus1        cell
        fmin          (1,1) double
        fmax          (1,1) double
        color
        options.Theme         string  = "default"
        options.Mirror        (1,1) logical = true
        options.gminus        cell   = {}
        options.Fill          (1,1) logical = false
        options.FaceAlpha     (1,1) double  = 0.3
        options.ShowBoundary  (1,1) logical = true
        options.BoundaryColor              = []    % [] → same as fill color
        options.FillTightness (1,1) double  = 0    % 0 = convex, 1 = exact
    end

    s = srg_style('Theme', options.Theme);
    c = resolve_color(color, s);

    % Resolve boundary color (defaults to fill color when not specified)
    if isempty(options.BoundaryColor)
        bc = c;
    else
        bc = resolve_color(options.BoundaryColor, s);
    end

    % Clamp FillTightness and snap to 0 / 1
    tightness = max(0, min(1, options.FillTightness));
    use_exact = (tightness >= 0.5);   % false → convex hull, true → polyshape

    % Frequency range indices
    if fmin == 0, initial = 1;
    else, [~, initial] = min(abs(fmin - rawdata.freq')); end
    if fmax == 0, final = length(rawdata.freq);
    else, [~, final] = min(abs(fmax - rawdata.freq')); end

    idx_range = initial:final;

    % Detect SISO: each slice is a single distinct value (srg_compute stores
    % field_of_values of a 1x1 matrix as a vector of identical scalars).
    is_siso = true;
    tol_det = 1e-10;
    for ii = idx_range
        cv = gplus1{ii};
        if isempty(cv), continue; end
        if length(cv) > 1 && max(abs(cv - cv(1))) > tol_det * max(abs(cv(1)), 1)
            is_siso = false;
            break;
        end
    end

    hold on

    % ------------------------------------------------------------------
    if is_siso

        if options.Fill
            warning('srg_plot_3d:FillNotApplicable', ...
                ['Fill = true is not applicable for SISO systems: each ' ...
                 'frequency slice is a single point, not a closed curve. ' ...
                 'The Fill option will be ignored.']);
        end

        % ---- SISO: collect scalars into 3D traces ----
        n = length(idx_range);
        gp_re  = zeros(n, 1);
        gp_im  = zeros(n, 1);
        freq_z = zeros(n, 1);

        for kk = 1:n
            ii = idx_range(kk);
            gp_re(kk)  = real(gplus1{ii}(1));
            gp_im(kk)  = imag(gplus1{ii}(1));
            freq_z(kk) = rawdata.freq(ii);
        end

        lw   = max(1.5, s.linewidth);
        mIdx = round(linspace(1, n, min(n, 25)));
        if n == 1
            plot3(gp_re, gp_im, freq_z, ...
                  'Color', c,   'LineStyle', 'none');
            if options.Mirror
                plot3(gp_re, -gp_im, freq_z, ...
                      'Color', c,   'LineStyle', 'none');
            end
        else
            % Line
            plot3(gp_re, gp_im, freq_z, ...
                  'Color', c, 'LineWidth', lw, 'LineStyle', '-');
            % Sparse markers for visibility (~25 evenly spaced)
            plot3(gp_re(mIdx), gp_im(mIdx), freq_z(mIdx), ...
                  'Color', c,   'LineStyle', 'none');
            if options.Mirror
                plot3(gp_re, -gp_im, freq_z, ...
                      'Color', c, 'LineWidth', lw, 'LineStyle', '-');
                plot3(gp_re(mIdx), -gp_im(mIdx), freq_z(mIdx), ...
                      'Color', c,   'LineStyle', 'none');
            end
        end

        % Also plot gminus trace if provided
        if ~isempty(options.gminus)
            gm_re = zeros(n, 1);
            gm_im = zeros(n, 1);
            for kk = 1:n
                ii = idx_range(kk);
                gm_re(kk) = real(options.gminus{ii}(1));
                gm_im(kk) = imag(options.gminus{ii}(1));
            end

            lw2   = max(1.5, s.linewidth);
            mIdx2 = round(linspace(1, n, min(n, 25)));
            if n == 1
                plot3(gm_re, gm_im, freq_z, ...
                      'Color', c,   'LineStyle', 'none');
                if options.Mirror
                    plot3(gm_re, -gm_im, freq_z, ...
                          'Color', c,   'LineStyle', 'none');
                end
            else
                plot3(gm_re, gm_im, freq_z, ...
                      'Color', c, 'LineWidth', lw2, 'LineStyle', '-');
                plot3(gm_re(mIdx2), gm_im(mIdx2), freq_z(mIdx2), ...
                      'Color', c,   'LineStyle', 'none');
                if options.Mirror
                    plot3(gm_re, -gm_im, freq_z, ...
                          'Color', c, 'LineWidth', lw2, 'LineStyle', '-');
                    plot3(gm_re(mIdx2), -gm_im(mIdx2), freq_z(mIdx2), ...
                          'Color', c,  'LineStyle', 'none');
                end
            end
        end

    % ------------------------------------------------------------------
    else
        % ---- MIMO: per-frequency curves (and optionally filled patches) ----

        for ii = idx_range
            freq_val = rawdata.freq(ii);
            curve = gplus1{ii};
            curve = curve(:);
            curve = curve(isfinite(curve));

            if isempty(curve), continue; end

            npts = length(curve);
            zz   = freq_val * ones(npts, 1);

            if npts == 1
                % Single point at this frequency: always use a marker
                plot3(real(curve), imag(curve), zz, 'o', ...
                      'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 10);
                if options.Mirror
                    plot3(real(curve), -imag(curve), zz, 'o', ...
                          'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 10);
                end

            else
                % ---- Filled patch ----
                if options.Fill
                    fill_freq_slice(real(curve), imag(curve), freq_val, ...
                                    c, options.FaceAlpha, use_exact);
                    if options.Mirror
                        fill_freq_slice(real(curve), -imag(curve), freq_val, ...
                                        c, options.FaceAlpha, use_exact);
                    end
                end

                % ---- Boundary wire ----
                % Draw when: (a) Fill is off, or (b) Fill is on and ShowBoundary is true.
                if ~options.Fill || options.ShowBoundary
                    wire_color = c;                  % default: fill/main color
                    if options.Fill
                        wire_color = bc;             % user-specified boundary color
                    end

                    plot3(real(curve), imag(curve), zz, ...
                          'Color', wire_color, 'LineWidth', s.linewidth);
                    if options.Mirror
                        plot3(real(curve), -imag(curve), zz, ...
                              'Color', wire_color, 'LineWidth', s.linewidth);
                    end
                end
            end
        end
    end

    % ------------------------------------------------------------------
    set(gca, 'ZScale', 'log');
    set(gca, 'FontSize', s.fontsize, 'FontName', s.fontname, ...
             'Color', s.bgcolor);
    grid(s.grid)
    box on
    xlabel('Re',        'FontSize', s.fontsize, 'FontName', s.fontname, ...
           'Interpreter', s.interpreter)
    ylabel('Im',        'FontSize', s.fontsize, 'FontName', s.fontname, ...
           'Interpreter', s.interpreter)
    zlabel('Frequency', 'FontSize', s.fontsize, 'FontName', s.fontname, ...
           'Interpreter', s.interpreter)
    view([-35, 25]);

end

%--------------------------------------------------------------------------
function fill_freq_slice(x, y, z_val, fc, fa, use_exact)
%FILL_FREQ_SLICE  Fill a single 2D slice at elevation z_val with fill3.
%
%   x, y     - Real coordinates of the ordered boundary curve (column vectors)
%   z_val    - Scalar frequency value (z elevation)
%   fc       - Face color (RGB triplet)
%   fa       - Face alpha in [0,1]
%   use_exact - false → convex hull; true → exact polyshape boundary

    if numel(x) < 3, return; end

    if ~use_exact
        % ---- Convex hull ----
        try
            k = convhull(x, y, 'Simplify', true);
            fill3(x(k), y(k), z_val * ones(size(k)), fc, ...
                  'FaceAlpha', fa, 'EdgeColor', 'none');
        catch
            % Degenerate slice (collinear points): silently skip
        end

    else
        % ---- Exact boundary via polyshape ----
        % The curve is already ordered, so we pass it directly.
        % 'Simplify', false preserves the shape without auto-convexification.
        try
            ps = polyshape(x, y, 'Simplify', false);
            if ps.NumRegions < 1, return; end
            [xb, yb] = boundary(ps);
            fill3(xb, yb, z_val * ones(size(xb)), fc, ...
                  'FaceAlpha', fa, 'EdgeColor', 'none');
        catch
            % Fallback to convex hull if polyshape construction fails
            try
                k = convhull(x, y, 'Simplify', true);
                fill3(x(k), y(k), z_val * ones(size(k)), fc, ...
                      'FaceAlpha', fa, 'EdgeColor', 'none');
            catch
            end
        end
    end

end

%--------------------------------------------------------------------------
function c = resolve_color(color, s)
%RESOLVE_COLOR  Always returns a plain 1x3 double RGB triplet.
%   Reads s.palette directly — avoids calling the s.color function handle,
%   which can behave unexpectedly inside nested functions on some versions.
    if isnumeric(color) && isscalar(color)
        n   = size(s.palette, 1);
        idx = mod(floor(color) - 1, n) + 1;
        c   = s.palette(idx, :);
    elseif isnumeric(color) && numel(color) == 3
        c = double(reshape(color, 1, 3));
    else
        c = s.palette(1, :);
    end
end