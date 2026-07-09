function srg_plot_3d_surface(rawdata, gplus, fmin, fmax, options)
%SRG_PLOT_3D_SURFACE  3D surface plot of SRG boundaries across frequencies.
%
%   SRG_PLOT_3D_SURFACE(RAWDATA, GPLUS, FMIN, FMAX) renders the SRG
%   as a continuous surface in 3D (Real, Imag, Frequency). Each frequency
%   slice is resampled to a common number of points, preserving the
%   original boundary traversal order from srg_compute.
%
%   Name-Value Arguments:
%       Theme       - Style theme (default: "default")
%       Color       - Color index (1-7) or RGB triplet (default: 1)
%       FaceAlpha   - Surface transparency (default: 0.6)
%       EdgeAlpha   - Edge transparency (default: 0.1)
%       Mirror      - Plot conjugate surface (default: true)
%       NumPoints   - Resampling points per curve (default: 100)
%       Lighting    - Enable lighting effects (default: true)
%       MinDiameter - Skip slices smaller than this (default: auto)
%
%   See also SRG_COMPUTE, SRG_PLOT_3D, SRG_PLOT_3D_COMPARE, SRG_HOMOTOPY

    arguments
        rawdata     table
        gplus       cell
        fmin        (1,1) double = 0
        fmax        (1,1) double = 0
        options.Theme       string = "default"
        options.Color       = 1
        options.FaceAlpha   (1,1) double = 0.6
        options.EdgeAlpha   (1,1) double = 0.1
        options.Mirror      (1,1) logical = true
        options.NumPoints   (1,1) double {mustBePositive, mustBeInteger} = 100
        options.Lighting    (1,1) logical = true
        options.MinDiameter (1,1) double = -1
    end

    s = srg_style('Theme', options.Theme);
    N = options.NumPoints;

    % Resolve color
    if isnumeric(options.Color) && isscalar(options.Color)
        face_color = s.color(options.Color);
    elseif isnumeric(options.Color) && length(options.Color) == 3
        face_color = options.Color;
    else
        face_color = s.color(1);
    end

    % Frequency range
    if fmin == 0, idx_start = 1;
    else, [~, idx_start] = min(abs(fmin - rawdata.freq')); end
    if fmax == 0, idx_end = length(rawdata.freq);
    else, [~, idx_end] = min(abs(fmax - rawdata.freq')); end

    freq_idx = idx_start:idx_end;
    n_freq = length(freq_idx);

    if n_freq < 2
        warning('srg_plot_3d_surface:tooFewFreqs', ...
                'Need at least 2 frequency points for a surface.');
        return;
    end

    % =====================================================================
    % Resample each slice preserving original point order
    % =====================================================================
    X = NaN(n_freq, N);
    Y = NaN(n_freq, N);
    Z = NaN(n_freq, N);
    diameters = zeros(n_freq, 1);
    valid_slice = false(n_freq, 1);

    for kk = 1:n_freq
        ii = freq_idx(kk);
        curve = gplus{ii};
        curve = curve(isfinite(curve));
        curve = curve(:);

        if length(curve) < 3
            continue;
        end

        % Diameter for filtering
        xr = real(curve); yr = imag(curve);
        diameters(kk) = sqrt((max(xr)-min(xr))^2 + (max(yr)-min(yr))^2);

        % Remove consecutive near-duplicates (preserve order)
        tol = max(1e-15, diameters(kk) * 1e-10);
        keep = [true; abs(diff(curve)) > tol];
        curve = curve(keep);

        if length(curve) < 3
            continue;
        end

        % Close the curve if not already closed
        if abs(curve(end) - curve(1)) > tol
            curve(end+1) = curve(1); %#ok<AGROW>
        end

        % Arc-length parameterization along the ORIGINAL order
        xc = real(curve);
        yc = imag(curve);
        ds = sqrt(diff(xc).^2 + diff(yc).^2);
        cumlen = [0; cumsum(ds)];

        % Remove zero-length segments
        good = [true; ds > tol];
        cumlen = cumlen(good);
        xc = xc(good);
        yc = yc(good);

        if length(cumlen) < 2 || cumlen(end) < tol
            continue;
        end

        % Resample to N equally spaced arc-length points
        target_s = linspace(0, cumlen(end), N);
        X(kk, :) = interp1(cumlen, xc, target_s, 'pchip');
        Y(kk, :) = interp1(cumlen, yc, target_s, 'pchip');
        Z(kk, :) = rawdata.freq(ii);
        valid_slice(kk) = true;
    end

    % =====================================================================
    % Filter degenerate slices
    % =====================================================================
    if options.MinDiameter < 0
        nonzero_d = diameters(valid_slice & diameters > 0);
        if ~isempty(nonzero_d)
            min_diam = 0.005 * median(nonzero_d);
        else
            min_diam = 0;
        end
    else
        min_diam = options.MinDiameter;
    end

    valid_slice = valid_slice & diameters >= min_diam;
    X = X(valid_slice, :);
    Y = Y(valid_slice, :);
    Z = Z(valid_slice, :);

    if size(X, 1) < 2
        warning('srg_plot_3d_surface:noValidSlices', ...
                'Too few valid frequency slices.');
        return;
    end

    % =====================================================================
    % Align starting point across slices for smooth stitching
    % =====================================================================
    % The original point order is consistent within each slice, but the
    % "starting point" on the closed curve can drift between frequencies.
    % Fix by circularly shifting each row to minimize distance to the
    % previous row's starting point.
    for kk = 2:size(X, 1)
        ref_x = X(kk-1, 1);
        ref_y = Y(kk-1, 1);
        dists = (X(kk,:) - ref_x).^2 + (Y(kk,:) - ref_y).^2;
        [~, best] = min(dists);
        if best > 1
            X(kk,:) = circshift(X(kk,:), -(best-1));
            Y(kk,:) = circshift(Y(kk,:), -(best-1));
        end
    end

    % =====================================================================
    % Render
    % =====================================================================
    hold on

    surf(X, Y, Z, ...
         'FaceColor', face_color, 'FaceAlpha', options.FaceAlpha, ...
         'EdgeColor', face_color, 'EdgeAlpha', options.EdgeAlpha);

    if options.Mirror
        surf(X, -Y, Z, ...
             'FaceColor', face_color, 'FaceAlpha', options.FaceAlpha, ...
             'EdgeColor', face_color, 'EdgeAlpha', options.EdgeAlpha);
    end

    set(gca, 'ZScale', 'log');

    if options.Lighting
        lighting gouraud
        material dull
        camlight('headlight');
        camlight('left');
    end

    set(gca, 'FontSize', s.fontsize, 'FontName', s.fontname, ...
             'Color', s.bgcolor);
    grid on
    box on
    xlabel('Re', 'FontSize', s.fontsize, 'FontName', s.fontname, ...
           'Interpreter', s.interpreter);
    ylabel('Im', 'FontSize', s.fontsize, 'FontName', s.fontname, ...
           'Interpreter', s.interpreter);
    zlabel('Frequency', 'FontSize', s.fontsize, 'FontName', s.fontname, ...
           'Interpreter', s.interpreter);

    view([35, 25]);

end
