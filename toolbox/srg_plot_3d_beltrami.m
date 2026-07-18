function srg_plot_3d_beltrami(rawdata, gplus1, fmin, fmax, color, options)
%SRG_PLOT_3D_BELTRAMI  3D wire plot in Beltrami-Klein plane across frequencies.
%
%   SRG_PLOT_3D_BELTRAMI(rawdata, gplus1, fmin, fmax, color)
%
%   Inputs:
%       rawdata - Table with .freq field
%       gplus1  - Cell array of BK-plane boundary data per frequency
%       fmin    - Minimum frequency (0 = first available)
%       fmax    - Maximum frequency (0 = last available)
%       color   - Color index (1-7) or RGB triplet
%
%   Name-Value Arguments:
%       Theme       - Style theme (default: "default")
%       ShowDisk    - Render the unit-circle as a translucent cylinder
%                     spanning the plotted frequency range (default: true)
%       DiskAlpha   - Disk surface transparency (default: 0.06)
%       DiskColor   - Disk color, index or RGB (default: [], neutral gray
%                     regardless of theme)
%       DiskNTheta  - Angular resolution of the cylinder (default: 180)
%
%   Behaviour:
%       SISO  — each cell contains one point; plotted as a continuous 3D
%               trajectory connecting points across frequencies.
%       MIMO  — each cell contains a curve; plotted as frequency slices.
%       Const — single-cell input (static matrix); curve replicated at
%               every frequency level.
%
%   See also SRG_COMPUTE, SRG_PLOT_3D, SRG_PLOT_3D_BELTRAMI_COMPARE
%
%   Example:
%       srg_plot_3d_beltrami(rawdata, W, 0, 0, 1)
%       srg_plot_3d_beltrami(rawdata, W, 0, 0, 1, ...
%           'DiskColor', [0.3 0.3 0.3], 'DiskAlpha', 0.1)
%       srg_plot_3d_beltrami(rawdata, W, 0, 0, 1, 'ShowDisk', false)

    arguments
        rawdata
        gplus1  cell
        fmin    (1,1) double
        fmax    (1,1) double
        color
        options.Theme string = "default"
        options.ShowDisk   (1,1) logical = true
        options.DiskAlpha  (1,1) double {mustBeInRange(options.DiskAlpha,0,1)} = 0.06
        options.DiskColor         = []
        options.DiskNTheta (1,1) double {mustBePositive, mustBeInteger} = 180
    end

    s = srg_style('Theme', options.Theme);
    c = resolve_color(color, s);

    if isempty(options.DiskColor)
        cd = [0.55 0.55 0.55];   % neutral gray regardless of theme
    else
        cd = resolve_color(options.DiskColor, s);
    end

    % --- frequency window ---------------------------------------------------
    if fmin == 0, initial = 1;
    else, [~, initial] = min(abs(fmin - rawdata.freq')); end
    if fmax == 0, final = length(rawdata.freq);
    else, [~, final] = min(abs(fmax - rawdata.freq')); end

    n_cells = numel(gplus1);

    hold on

    % ========================================================================
    % Case 1: CONSTANT — caller passed a single-cell {W_matrix}
    % ========================================================================
    if n_cells == 1
        crv = gplus1{1};
        npts = length(crv);
        for ii = initial:final
            fv = rawdata.freq(ii);
            plot3(real(crv), imag(crv), ones(npts,1)*fv, ...
                  'Color', c, 'LineWidth', s.linewidth);
        end

    % ========================================================================
    % Case 2: SISO — every cell is a single unique point (scalar TF / 1×1)
    % ========================================================================
    elseif is_siso_data(gplus1, initial, final)
        n_pts = final - initial + 1;
        traj  = zeros(n_pts, 1);
        freqs = zeros(n_pts, 1);
        for ii = initial:final
            crv = gplus1{ii}(isfinite(gplus1{ii}));
            traj(ii-initial+1)  = crv(1);          % one representative point
            freqs(ii-initial+1) = rawdata.freq(ii);
        end
        plot3(real(traj), imag(traj), freqs, ...
              '-', 'Color', c, 'LineWidth', s.linewidth + 0.5);
        plot3(real(traj), imag(traj), freqs, ...
              '.', 'Color', c, 'MarkerSize', 4);     % frequency markers

    % ========================================================================
    % Case 3: MIMO — each cell is a curve; plot per-frequency slice
    % ========================================================================
    else
        for ii = initial:final
            fv   = rawdata.freq(ii);
            crv  = gplus1{ii};
            npts = length(crv);

            % Degenerate slice: all points nearly identical → dot marker
            if is_single_point(crv)
                plot3(real(crv(1)), imag(crv(1)), fv, ...
                      '.', 'Color', c, 'MarkerSize', 6);
            else
                plot3(real(crv), imag(crv), ones(npts,1)*fv, ...
                      'Color', c, 'LineWidth', s.linewidth);
            end
        end
    end

    % ========================================================================
    % Unit-disk cylinder boundary
    % ========================================================================
    if options.ShowDisk
        f_lo = min(rawdata.freq(initial), rawdata.freq(final));
        f_hi = max(rawdata.freq(initial), rawdata.freq(final));
        draw_disk_cylinder(f_lo, f_hi, cd, options.DiskAlpha, options.DiskNTheta);
    end

    set(gca, 'ZScale', 'log');
    set(gca, 'FontSize', s.fontsize, 'FontName', s.fontname, ...
             'Color', s.bgcolor);
    grid(s.grid)
    box on
    xlabel('Real',      'FontSize', s.fontsize, 'FontName', s.fontname, ...
                        'Interpreter', s.interpreter)
    ylabel('Imag',      'FontSize', s.fontsize, 'FontName', s.fontname, ...
                        'Interpreter', s.interpreter)
    zlabel('Frequency', 'FontSize', s.fontsize, 'FontName', s.fontname, ...
                        'Interpreter', s.interpreter)
end


% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function tf = is_siso_data(gp, i0, i1)
%IS_SISO_DATA  True when every cell in the window contains one unique point.
    tf = true;
    for ii = i0:i1
        if ~is_single_point(gp{ii})
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

function draw_disk_cylinder(f_lo, f_hi, cd, disk_alpha, ntheta)
%DRAW_DISK_CYLINDER  Render the unit-circle as a translucent 3D cylinder.
%   Uses a log-spaced set of z levels so it looks smooth on a log z-axis.

    % Guard against a non-positive or degenerate frequency window (e.g.
    % a single-frequency constant-matrix plot, or a linear sweep that
    % starts at 0 Hz): log10 of a non-positive value is undefined, and a
    % single point would otherwise draw a zero-height cylinder anyway.
    if f_lo <= 0
        f_lo = max(f_hi * 1e-6, eps);
    end
    if f_hi <= f_lo
        f_hi = f_lo * (1 + 1e-6);
    end

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

    % Subtle top and bottom caps (circles at the frequency limits)
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
