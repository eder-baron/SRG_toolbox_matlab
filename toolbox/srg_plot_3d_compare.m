function srg_plot_3d_compare(rawdata1, gplus1, rawdata2, gplus2, fmin, fmax, options)
%SRG_PLOT_3D_COMPARE  Overlay two SRG surfaces in 3D with intersection highlighting.
%
%   SRG_PLOT_3D_COMPARE(RAWDATA1, GPLUS1, RAWDATA2, GPLUS2, FMIN, FMAX)
%   renders two SRG surfaces in 3D and highlights frequency slices where
%   the two SRGs overlap.
%
%   Name-Value Arguments:
%       Theme          - Style theme (default: "default")
%       Color1         - Color for SRG 1 (index or RGB, default: 1)
%       Color2         - Color for SRG 2 (index or RGB, default: 5)
%       IntersectColor - Color for intersection (default: CB red)
%       FaceAlpha      - Surface face transparency (default: 0.45)
%       EdgeAlpha      - Surface edge-line opacity  (default: 0.50)
%       IntersectAlpha - Intersection patch opacity (default: 0.9)
%       Mirror         - Plot conjugate surfaces (default: true)
%       NumPoints      - Resampling points per curve (default: 100)
%       Lighting       - Enable lighting (default: true)
%       MinDiameter    - Skip tiny slices (default: auto)
%
%   Tip: if surfaces still look too transparent, raise FaceAlpha toward 1
%   and/or set EdgeAlpha to 0.6–0.8 so the wireframe grid reads clearly.
%
%   See also SRG_PLOT_3D_SURFACE, SRG_COMPUTE, SRG_STABILITY_MARGIN, SRG_PLOT_WITNESS_3D

    arguments
        rawdata1    table
        gplus1      cell
        rawdata2    table
        gplus2      cell
        fmin        (1,1) double = 0
        fmax        (1,1) double = 0
        options.Theme          string  = "default"
        options.Color1                 = 1
        options.Color2                 = 5
        options.IntersectColor         = []
        options.FaceAlpha      (1,1) double = 0.45
        options.EdgeAlpha      (1,1) double = 0.50   % <-- was hardcoded 0.08
        options.IntersectAlpha (1,1) double = 0.9
        options.Mirror         (1,1) logical = true
        options.NumPoints      (1,1) double {mustBePositive, mustBeInteger} = 100
        options.Lighting       (1,1) logical = true
        options.MinDiameter    (1,1) double = -1
    end

    s  = srg_style('Theme', options.Theme);
    N  = options.NumPoints;
    ea = options.EdgeAlpha;

    % Resolve colors (pass s so theme is respected)
    c1 = resolve_color(options.Color1, s);
    c2 = resolve_color(options.Color2, s);
    if isempty(options.IntersectColor)
        ci = s.cb.red;
    else
        ci = resolve_color(options.IntersectColor, s);
    end

    % =====================================================================
    % SISO short-circuit
    % =====================================================================
    if is_all_constant(gplus1) && is_all_constant(gplus2)
        srg_plot_3d(rawdata1, gplus1, fmin, fmax, options.Color1, ...
                    'Theme', options.Theme, 'Mirror', options.Mirror);
        srg_plot_3d(rawdata2, gplus2, fmin, fmax, options.Color2, ...
                    'Theme', options.Theme, 'Mirror', options.Mirror);
        view([-35, 25]);
        grid on;
        return;
    end

    % =====================================================================
    % Find common frequency grid
    % =====================================================================
    freq1 = rawdata1.freq;
    freq2 = rawdata2.freq;

    if fmin == 0, f_lo = max(min(freq1), min(freq2));
    else,         f_lo = fmin; end
    if fmax == 0, f_hi = min(max(freq1), max(freq2));
    else,         f_hi = fmax; end

    mask1  = freq1 >= f_lo & freq1 <= f_hi;
    idx1   = find(mask1);
    n_freq = length(idx1);

    if n_freq < 2
        warning('srg_plot_3d_compare:tooFewFreqs', ...
                'Need at least 2 common frequency points.');
        return;
    end

    idx2 = zeros(n_freq, 1);
    for kk = 1:n_freq
        [~, idx2(kk)] = min(abs(freq2 - freq1(idx1(kk))));
    end

    % =====================================================================
    % Resample both systems and detect intersections
    % =====================================================================
    X1 = NaN(n_freq, N); Y1 = NaN(n_freq, N); Z1 = NaN(n_freq, N);
    X2 = NaN(n_freq, N); Y2 = NaN(n_freq, N); Z2 = NaN(n_freq, N);
    valid1 = false(n_freq, 1);
    valid2 = false(n_freq, 1);

    intersect_patches = cell(n_freq, 1);
    intersect_freqs   = NaN(n_freq, 1);

    for kk = 1:n_freq
        f = freq1(idx1(kk));

        [xr1, yr1, ok1] = resample_slice(gplus1{idx1(kk)}, N);
        if ok1
            X1(kk,:) = xr1; Y1(kk,:) = yr1; Z1(kk,:) = f;
            valid1(kk) = true;
        end

        [xr2, yr2, ok2] = resample_slice(gplus2{idx2(kk)}, N);
        if ok2
            X2(kk,:) = xr2; Y2(kk,:) = yr2; Z2(kk,:) = f;
            valid2(kk) = true;
        end

        if ok1 && ok2
            ps1 = make_polyshape(gplus1{idx1(kk)});
            ps2 = make_polyshape(gplus2{idx2(kk)});
            if ps1.NumRegions > 0 && ps2.NumRegions > 0
                ps_int = intersect(ps1, ps2);
                if ps_int.NumRegions > 0 && area(ps_int) > 0
                    intersect_patches{kk} = ps_int;
                    intersect_freqs(kk)   = f;
                end
            end
        end
    end

    % Auto min diameter filter
    diam1 = compute_diameters(X1, Y1, valid1);
    diam2 = compute_diameters(X2, Y2, valid2);

    if options.MinDiameter < 0
        all_d = [diam1(valid1 & diam1>0); diam2(valid2 & diam2>0)];
        min_diam = 0;
        if ~isempty(all_d)
            min_diam = 0.005 * median(all_d);
        end
    else
        min_diam = options.MinDiameter;
    end

    valid1 = valid1 & diam1 >= min_diam;
    valid2 = valid2 & diam2 >= min_diam;

    % =====================================================================
    % Align starting points across slices
    % =====================================================================
    [X1, Y1] = align_slices(X1, Y1, valid1);
    [X2, Y2] = align_slices(X2, Y2, valid2);

    % =====================================================================
    % Render
    % =====================================================================
    hold on

    % Tell MATLAB to sort transparent faces by depth — fixes the
    % additive-blending wash-out that occurs with multiple stacked surfaces
    set(gca, 'SortMethod', 'depth');

    % Surface 1
    if sum(valid1) >= 2
        surf(X1(valid1,:), Y1(valid1,:), Z1(valid1,:), ...
             'FaceColor', c1, 'FaceAlpha', options.FaceAlpha, ...
             'EdgeColor', c1, 'EdgeAlpha', ea);
        if options.Mirror
            surf(X1(valid1,:), -Y1(valid1,:), Z1(valid1,:), ...
                 'FaceColor', c1, 'FaceAlpha', options.FaceAlpha, ...
                 'EdgeColor', c1, 'EdgeAlpha', ea);
        end
    end

    % Surface 2
    if sum(valid2) >= 2
        surf(X2(valid2,:), Y2(valid2,:), Z2(valid2,:), ...
             'FaceColor', c2, 'FaceAlpha', options.FaceAlpha, ...
             'EdgeColor', c2, 'EdgeAlpha', ea);
        if options.Mirror
            surf(X2(valid2,:), -Y2(valid2,:), Z2(valid2,:), ...
                 'FaceColor', c2, 'FaceAlpha', options.FaceAlpha, ...
                 'EdgeColor', c2, 'EdgeAlpha', ea);
        end
    end

    % Intersection patches
    for kk = 1:n_freq
        if ~isempty(intersect_patches{kk})
            ps_int = intersect_patches{kk};
            f      = intersect_freqs(kk);
            [xv, yv] = boundary(ps_int);

            patch('XData', xv, 'YData', yv, ...
                  'ZData', f * ones(size(xv)), ...
                  'FaceColor', ci, 'FaceAlpha', options.IntersectAlpha, ...
                  'EdgeColor', ci, 'EdgeAlpha', 0.6, 'LineWidth', 1.2);

            if options.Mirror
                patch('XData', xv, 'YData', -yv, ...
                      'ZData', f * ones(size(xv)), ...
                      'FaceColor', ci, 'FaceAlpha', options.IntersectAlpha, ...
                      'EdgeColor', ci, 'EdgeAlpha', 0.6, 'LineWidth', 1.2);
            end
        end
    end

    % =====================================================================
    % Axes styling
    % =====================================================================
    set(gca, 'ZScale', 'log');

    if options.Lighting
        lighting gouraud
        material dull
        camlight('headlight');
        camlight('left');
    end

    set(gca, 'FontSize', s.fontsize, 'FontName', s.fontname, ...
             'Color', s.bgcolor);
    grid on; box on;
    xlabel('Re', 'FontSize', s.fontsize, 'FontName', s.fontname, ...
           'Interpreter', s.interpreter);
    ylabel('Im', 'FontSize', s.fontsize, 'FontName', s.fontname, ...
           'Interpreter', s.interpreter);
    zlabel('Frequency', 'FontSize', s.fontsize, 'FontName', s.fontname, ...
           'Interpreter', s.interpreter);
    view([-35, 25]);

end


% =========================================================================
% Local helper functions
% =========================================================================

function c = resolve_color(color, s)
%RESOLVE_COLOR  Convert color index or RGB to an RGB triplet using style s.
    if isnumeric(color) && isscalar(color)
        c = s.color(color);             % use the theme palette
    elseif isnumeric(color) && numel(color) == 3
        c = color(:)';
    else
        c = s.color(1);
    end
end


function [xr, yr, ok] = resample_slice(curve, N)
%RESAMPLE_SLICE  Resample preserving original boundary order.
    ok = false;
    xr = NaN(1, N);
    yr = NaN(1, N);

    curve = curve(isfinite(curve));
    curve = curve(:);

    if length(curve) < 3, return; end

    tol  = max(1e-15, max(abs(curve)) * 1e-10);
    keep = [true; abs(diff(curve)) > tol];
    curve = curve(keep);

    if length(curve) < 3, return; end

    if abs(curve(end) - curve(1)) > tol
        curve(end+1) = curve(1);
    end

    xc = real(curve);
    yc = imag(curve);
    ds = sqrt(diff(xc).^2 + diff(yc).^2);
    cumlen = [0; cumsum(ds)];

    good = [true; ds > tol];
    cumlen = cumlen(good);
    xc = xc(good);
    yc = yc(good);

    if length(cumlen) < 2 || cumlen(end) < tol, return; end

    target_s = linspace(0, cumlen(end), N);
    xr = interp1(cumlen, xc, target_s, 'pchip');
    yr = interp1(cumlen, yc, target_s, 'pchip');
    ok = true;
end


function ps = make_polyshape(curve)
%MAKE_POLYSHAPE  Convert complex curve to polyshape preserving original order.
    curve = curve(isfinite(curve));
    curve = curve(:);

    if length(curve) < 3, ps = polyshape(); return; end

    tol  = max(1e-15, max(abs(curve)) * 1e-10);
    keep = [true; abs(diff(curve)) > tol];
    curve = curve(keep);

    if length(curve) < 3, ps = polyshape(); return; end

    ps = polyshape(real(curve), imag(curve), 'Simplify', false);
end


function [X, Y] = align_slices(X, Y, valid)
%ALIGN_SLICES  Circularly shift rows to align starting points.
    valid_idx = find(valid);
    for kk = 2:length(valid_idx)
        row  = valid_idx(kk);
        prev = valid_idx(kk-1);
        dists = (X(row,:) - X(prev,1)).^2 + (Y(row,:) - Y(prev,1)).^2;
        [~, best] = min(dists);
        if best > 1
            X(row,:) = circshift(X(row,:), -(best-1));
            Y(row,:) = circshift(Y(row,:), -(best-1));
        end
    end
end


function diams = compute_diameters(X, Y, valid)
    n     = size(X, 1);
    diams = zeros(n, 1);
    for kk = 1:n
        if valid(kk)
            diams(kk) = sqrt((max(X(kk,:))-min(X(kk,:)))^2 + ...
                             (max(Y(kk,:))-min(Y(kk,:)))^2);
        end
    end
end


function tf = is_all_constant(gp_cells)
%IS_ALL_CONSTANT  True iff every cell entry is a single distinct value.
    tf  = true;
    tol = 1e-10;
    for ii = 1:numel(gp_cells)
        c = gp_cells{ii};
        if isempty(c), continue; end
        if length(c) > 1 && max(abs(c - c(1))) > tol * max(abs(c(1)), 1)
            tf = false;
            return;
        end
    end
end