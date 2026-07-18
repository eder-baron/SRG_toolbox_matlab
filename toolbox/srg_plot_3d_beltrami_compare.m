function srg_plot_3d_beltrami_compare(rawdata1, W1, rawdata2, W2, fmin, fmax, options)
%SRG_PLOT_3D_BELTRAMI_COMPARE  Overlay two SRG surfaces in the Beltrami-Klein
%   disk across frequency, with optional intersection highlighting and unit
%   disk boundary.
%
%   SRG_PLOT_3D_BELTRAMI_COMPARE(RD1, W1, RD2, W2, FMIN, FMAX)
%
%   Inputs:
%       rawdata1/2  - Tables with .freq field (from SRG_COMPUTE or SRG)
%       W1/W2       - Cell arrays of BK-plane boundary curves per frequency
%       fmin        - Minimum frequency to render  (0 = first available)
%       fmax        - Maximum frequency to render  (0 = last available)
%
%   Name-Value Arguments:
%       Theme          - Style theme                      (default: "default")
%       Color1         - Color for W1 (index or RGB)      (default: 1)
%       Color2         - Color for W2 (index or RGB)      (default: 5)
%       IntersectColor - Color for overlap patches        (default: CB red)
%       FaceAlpha      - Surface transparency             (default: 0.45)
%       IntersectAlpha - Intersection patch opacity       (default: 0.9)
%       Lighting       - Enable Gouraud lighting          (default: true)
%       NumPoints      - Resampling points per curve      (default: 100)
%       MinDiameter    - Skip slices thinner than this    (default: auto)
%       ShowDisk       - Render the unit-circle cylinder  (default: true)
%       DiskAlpha      - Transparency of disk surface     (default: 0.06)
%       DiskColor      - Color of disk surface/edge       (default: gray)
%       DiskNTheta     - Angular resolution of cylinder   (default: 180)
%
%   See also SRG_PLOT_3D_BELTRAMI, SRG_PLOT_3D_COMPARE, SRG_COMPUTE

    arguments
        rawdata1   table
        W1         cell
        rawdata2   table
        W2         cell
        fmin       (1,1) double = 0
        fmax       (1,1) double = 0
        options.Theme          string  = "default"
        options.Color1                 = 1
        options.Color2                 = 5
        options.IntersectColor         = []
        options.FaceAlpha      (1,1) double = 0.45
        options.IntersectAlpha (1,1) double = 0.9
        options.Lighting       (1,1) logical = true
        options.NumPoints      (1,1) double {mustBePositive,mustBeInteger} = 100
        options.MinDiameter    (1,1) double = -1
        options.ShowDisk       (1,1) logical = true
        options.DiskAlpha      (1,1) double = 0.06
        options.DiskColor              = []
        options.DiskNTheta     (1,1) double {mustBePositive,mustBeInteger} = 180
    end

    s  = srg_style('Theme', options.Theme);
    N  = options.NumPoints;

    c1 = resolve_color(options.Color1, s);
    c2 = resolve_color(options.Color2, s);

    if isempty(options.IntersectColor)
        ci = s.cb.red;
    else
        ci = resolve_color(options.IntersectColor, s);
    end

    if isempty(options.DiskColor)
        cd = [0.55 0.55 0.55];   % neutral gray regardless of theme
    else
        cd = resolve_color(options.DiskColor, s);
    end

    % ------------------------------------------------------------------
    % Frequency range
    % ------------------------------------------------------------------
    freq1 = rawdata1.freq;
    freq2 = rawdata2.freq;

    if fmin == 0, f_lo = max(min(freq1), min(freq2));
    else,         f_lo = fmin; end
    if fmax == 0, f_hi = min(max(freq1), max(freq2));
    else,         f_hi = fmax; end

    mask1  = freq1 >= f_lo & freq1 <= f_hi;
    idx1   = find(mask1);
    n_freq = numel(idx1);

    if n_freq < 2
        warning('srg_plot_3d_beltrami_compare:tooFewFreqs', ...
                'Need at least 2 common frequency points.');
        return;
    end

    idx2 = zeros(n_freq, 1);
    for kk = 1:n_freq
        [~, idx2(kk)] = min(abs(freq2 - freq1(idx1(kk))));
    end

    % ------------------------------------------------------------------
    % SISO detection (per input — a nominal MIMO system can still collapse
    % to a single point per frequency, e.g. scalar * I)
    % ------------------------------------------------------------------
    siso1 = is_all_constant(W1);
    siso2 = is_all_constant(W2);

    if siso1 && siso2
        srg_plot_3d_beltrami(rawdata1, W1, fmin, fmax, options.Color1, ...
                             'Theme', options.Theme);
        srg_plot_3d_beltrami(rawdata2, W2, fmin, fmax, options.Color2, ...
                             'Theme', options.Theme);
        view([-35, 25]);
        grid(s.grid)
        if options.ShowDisk
            draw_disk_cylinder(f_lo, f_hi, cd, options.DiskAlpha, options.DiskNTheta);
        end
        return;
    end

    % ------------------------------------------------------------------
    % Mixed case: one side is a genuine MIMO surface, the other collapses
    % to a single point per frequency. A polygon/polygon intersection
    % against a degenerate curve is meaningless, so instead render the
    % degenerate side as a trajectory and check, at each frequency,
    % whether that point lies inside the MIMO SRG (point-in-polygon).
    % ------------------------------------------------------------------
    if siso1 ~= siso2
        if siso1
            pts_siso = siso_trajectory(W1, idx1);
            W_mimo   = W2;  idx_mimo = idx2;
            c_siso   = c1;  c_mimo   = c2;
        else
            pts_siso = siso_trajectory(W2, idx2);
            W_mimo   = W1;  idx_mimo = idx1;
            c_siso   = c2;  c_mimo   = c1;
        end

        freqs = freq1(idx1);

        Xm = NaN(n_freq, N);  Ym = NaN(n_freq, N);  Zm = NaN(n_freq, N);
        valid_m = false(n_freq, 1);
        inside  = false(n_freq, 1);

        for kk = 1:n_freq
            [xr, yr, ok] = resample_slice(W_mimo{idx_mimo(kk)}, N);
            if ok
                Xm(kk,:) = xr;  Ym(kk,:) = yr;  Zm(kk,:) = freqs(kk);
                valid_m(kk) = true;
            end

            if isfinite(pts_siso(kk))
                ps = make_polyshape(W_mimo{idx_mimo(kk)});
                if ps.NumRegions > 0
                    inside(kk) = isinterior(ps, real(pts_siso(kk)), imag(pts_siso(kk)));
                end
            end
        end

        diam_m = compute_diameters(Xm, Ym, valid_m);
        if options.MinDiameter < 0
            d = diam_m(valid_m & diam_m > 0);
            min_diam = 0;
            if ~isempty(d), min_diam = 0.005 * median(d); end
        else
            min_diam = options.MinDiameter;
        end
        valid_m = valid_m & diam_m >= min_diam;

        [Xm, Ym] = align_slices(Xm, Ym, valid_m);

        hold on

        if sum(valid_m) >= 2
            surf(Xm(valid_m,:), Ym(valid_m,:), Zm(valid_m,:), ...
                 'FaceColor', c_mimo, 'FaceAlpha', options.FaceAlpha, ...
                 'EdgeColor', c_mimo, 'EdgeAlpha', 0.08);
        end

        ok_pts = isfinite(pts_siso);
        plot3(real(pts_siso(ok_pts)), imag(pts_siso(ok_pts)), freqs(ok_pts), ...
              '-', 'Color', c_siso, 'LineWidth', s.linewidth);

        if any(inside)
            plot3(real(pts_siso(inside)), imag(pts_siso(inside)), freqs(inside), ...
                  'o', 'MarkerFaceColor', ci, 'MarkerEdgeColor', ci, 'MarkerSize', 5);
            warning('srg_plot_3d_beltrami_compare:trajectoryInsideSurface', ...
                'SISO trajectory enters the MIMO SRG at %d of %d frequencies (%.3g-%.3g Hz).', ...
                nnz(inside), n_freq, min(freqs(inside)), max(freqs(inside)));
        end

        if options.ShowDisk
            draw_disk_cylinder(f_lo, f_hi, cd, options.DiskAlpha, options.DiskNTheta);
        end

        apply_axes_styling(s, options);
        return;
    end

    % ------------------------------------------------------------------
    % Resample and detect intersections
    % ------------------------------------------------------------------
    X1 = NaN(n_freq, N);  Y1 = NaN(n_freq, N);  Z1 = NaN(n_freq, N);
    X2 = NaN(n_freq, N);  Y2 = NaN(n_freq, N);  Z2 = NaN(n_freq, N);
    valid1 = false(n_freq, 1);
    valid2 = false(n_freq, 1);

    intersect_patches = cell(n_freq, 1);
    intersect_freqs   = NaN(n_freq, 1);

    for kk = 1:n_freq
        f = freq1(idx1(kk));

        [xr1, yr1, ok1] = resample_slice(W1{idx1(kk)}, N);
        if ok1
            X1(kk,:) = xr1;  Y1(kk,:) = yr1;  Z1(kk,:) = f;
            valid1(kk) = true;
        end

        [xr2, yr2, ok2] = resample_slice(W2{idx2(kk)}, N);
        if ok2
            X2(kk,:) = xr2;  Y2(kk,:) = yr2;  Z2(kk,:) = f;
            valid2(kk) = true;
        end

        if ok1 && ok2
            ps1 = make_polyshape(W1{idx1(kk)});
            ps2 = make_polyshape(W2{idx2(kk)});
            if ps1.NumRegions > 0 && ps2.NumRegions > 0
                ps_int = intersect(ps1, ps2);
                if ps_int.NumRegions > 0 && area(ps_int) > 0
                    intersect_patches{kk} = ps_int;
                    intersect_freqs(kk)   = f;
                end
            end
        end
    end

    % Auto min-diameter filter
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

    % ------------------------------------------------------------------
    % Align starting points across slices
    % ------------------------------------------------------------------
    [X1, Y1] = align_slices(X1, Y1, valid1);
    [X2, Y2] = align_slices(X2, Y2, valid2);

    % ------------------------------------------------------------------
    % Render
    % ------------------------------------------------------------------
    hold on

    % Surface 1
    if sum(valid1) >= 2
        surf(X1(valid1,:), Y1(valid1,:), Z1(valid1,:), ...
             'FaceColor', c1, 'FaceAlpha', options.FaceAlpha, ...
             'EdgeColor', c1, 'EdgeAlpha', 0.08);
    end

    % Surface 2
    if sum(valid2) >= 2
        surf(X2(valid2,:), Y2(valid2,:), Z2(valid2,:), ...
             'FaceColor', c2, 'FaceAlpha', options.FaceAlpha, ...
             'EdgeColor', c2, 'EdgeAlpha', 0.08);
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
        end
    end

    % Unit-disk cylinder boundary
    if options.ShowDisk
        draw_disk_cylinder(f_lo, f_hi, cd, options.DiskAlpha, options.DiskNTheta);
    end

    % ------------------------------------------------------------------
    % Axes styling
    % ------------------------------------------------------------------
    apply_axes_styling(s, options);

end


% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function draw_disk_cylinder(f_lo, f_hi, cd, disk_alpha, ntheta)
%DRAW_DISK_CYLINDER  Render the unit-circle as a translucent 3D cylinder.
%   Uses a log-spaced set of z levels so it looks smooth on a log z-axis.

    n_z    = 60;
    theta  = linspace(0, 2*pi, ntheta);
    z_vals = logspace(log10(f_lo), log10(f_hi), n_z);

    Xd = cos(theta);      % 1 x ntheta
    Yd = sin(theta);      % 1 x ntheta

    % Replicate to n_z x ntheta grids
    Xd = repmat(Xd, n_z, 1);
    Yd = repmat(Yd, n_z, 1);
    Zd = repmat(z_vals(:), 1, ntheta);

    surf(Xd, Yd, Zd, ...
         'FaceColor', cd, 'FaceAlpha', disk_alpha, ...
         'EdgeColor', cd, 'EdgeAlpha', disk_alpha * 1.5, ...
         'LineStyle', '-', 'LineWidth', 0.3, ...
         'HandleVisibility', 'off');

    % Subtle top and bottom caps (circles at freq limits)
    for fz = [f_lo, f_hi]
        fill3(cos(theta), sin(theta), fz*ones(1,ntheta), cd, ...
              'FaceAlpha', disk_alpha * 0.8, ...
              'EdgeColor', cd, 'EdgeAlpha', disk_alpha * 2, ...
              'LineWidth', 0.5, 'HandleVisibility', 'off');
    end
end


function c = resolve_color(color, s)
    if isnumeric(color) && isscalar(color)
        c = s.color(color);
    elseif isnumeric(color) && numel(color) == 3
        c = color(:)';
    else
        c = s.color(1);
    end
end


function [xr, yr, ok] = resample_slice(curve, N)
    ok = false;
    xr = NaN(1, N);
    yr = NaN(1, N);

    curve = curve(isfinite(curve));
    curve = curve(:);

    if length(curve) < 3
        return;
    end

    tol  = max(1e-15, max(abs(curve)) * 1e-10);
    keep = [true; abs(diff(curve)) > tol];
    curve = curve(keep);

    if length(curve) < 3
        return;
    end

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

    if length(cumlen) < 2 || cumlen(end) < tol
        return;
    end

    target_s = linspace(0, cumlen(end), N);
    xr = interp1(cumlen, xc, target_s, 'pchip');
    yr = interp1(cumlen, yc, target_s, 'pchip');
    ok = true;
end


function ps = make_polyshape(curve)
    curve = curve(isfinite(curve));
    curve = curve(:);

    if length(curve) < 3
        ps = polyshape();
        return;
    end

    tol  = max(1e-15, max(abs(curve)) * 1e-10);
    keep = [true; abs(diff(curve)) > tol];
    curve = curve(keep);

    if length(curve) < 3
        ps = polyshape();
        return;
    end

    ps = polyshape(real(curve), imag(curve), 'Simplify', false);
end


function [X, Y] = align_slices(X, Y, valid)
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


function apply_axes_styling(s, options)
    set(gca, 'ZScale', 'log');

    if options.Lighting
        lighting gouraud
        material dull
        camlight('headlight');
        camlight('left');
    end

    set(gca, 'FontSize', s.fontsize, 'FontName', s.fontname, ...
             'Color', s.bgcolor);
    grid on;  box on;
    xlim([-1.05, 1.05]);
    ylim([-1.05, 1.05]);
    xlabel('Re', 'FontSize', s.fontsize, 'FontName', s.fontname, ...
           'Interpreter', s.interpreter);
    ylabel('Im', 'FontSize', s.fontsize, 'FontName', s.fontname, ...
           'Interpreter', s.interpreter);
    zlabel('Frequency', 'FontSize', s.fontsize, 'FontName', s.fontname, ...
           'Interpreter', s.interpreter);
    view([-35, 25]);
end


function pts = siso_trajectory(W, idx)
%SISO_TRAJECTORY  One representative complex point per frequency for a
%   cell array that has collapsed to a single point at every slice.
    n   = numel(idx);
    pts = NaN(n, 1);
    for kk = 1:n
        c = W{idx(kk)};
        c = c(isfinite(c));
        if ~isempty(c)
            pts(kk) = c(1);
        end
    end
end


function tf = is_all_constant(gp)
    tf  = true;
    tol = 1e-10;
    for ii = 1:numel(gp)
        c = gp{ii};
        if isempty(c), continue; end
        if length(c) > 1 && max(abs(c - c(1))) > tol * max(abs(c(1)), 1)
            tf = false;
            return;
        end
    end
end
