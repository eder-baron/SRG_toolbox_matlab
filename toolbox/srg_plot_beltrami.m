function srg_plot_beltrami(TransferFcn, W, A, color, estpoints, options)
%SRG_PLOT_BELTRAMI  Plot the SRG in the Beltrami-Klein disk.
%
%   SRG_PLOT_BELTRAMI(TransferFcn, W, A, color, estpoints)
%
%   Name-Value Arguments:
%       Theme     - Style theme (default: "default")
%       FaceAlpha - Fill transparency (default: 0.5)
%       Fill      - Fill the interior of each W slice (default: false)
%
%   See also SRG_COMPUTE, SRG_PLOT_GAUSS, SRG_PLOT_COMPARE

    arguments
        TransferFcn
        W           cell
        A           cell
        color
        estpoints   (1,1) double
        options.Theme     string = "default"
        options.FaceAlpha (1,1) double = 0.5
        options.Fill      (1,1) logical = false
    end

    s = srg_style('Theme', options.Theme);
    c = resolve_color(color, s);

    % Unit circle
    angles = linspace(0, 2*pi, 720);
    xc = cos(angles);
    yc = sin(angles);

    hold on

    % Detect SISO: the numerical range of a 1x1 matrix is a single point,
    % so W{i} is a vector of identical values. Treating that as a region
    % produces an invisible plot — instead, trace the scalar BK image
    % across frequency as a curve inside the disk.
    is_siso = ~isa(TransferFcn, 'double') && ...
              isequal(size(TransferFcn), [1 1]);
    if ~is_siso && ~isa(TransferFcn, 'double')
        is_siso = true;
        tol = 1e-10;
        for ii = 1:min(estpoints, length(W))
            wi = W{ii};
            if isempty(wi), continue; end
            if length(wi) > 1 && max(abs(wi - wi(1))) > tol * max(abs(wi(1)), 1)
                is_siso = false; break;
            end
        end
    end

    if isa(TransferFcn, 'double')
        if options.Fill
            fill_bk_slice(W{1}, c, options.FaceAlpha);
        end
        plot(W{1}, 'LineWidth', s.linewidth, 'Color', c)
       % plot(A{1}, 'x', 'LineWidth', s.linewidth, 'Color', c)
    elseif is_siso
        % --- SISO: BK trajectory of the scalar across frequency ---
        n = min(estpoints, length(W));
        W_tr = zeros(n, 1);
        for ii = 1:n
            W_tr(ii) = W{ii}(1);
        end
        mIdx = round(linspace(1, n, min(n, 25)));
        lw   = max(1.5, s.linewidth);
        plot(real(W_tr), imag(W_tr), '-', ...
             'Color', c, 'LineWidth', lw, ...
             'Marker', 'o', 'MarkerFaceColor', c, ...
             'MarkerEdgeColor', 'none', 'MarkerSize', 4, ...
             'MarkerIndices', mIdx);
    else
        for ii = 1:estpoints-1
            if options.Fill
                fill_bk_slice(W{ii}, c, options.FaceAlpha);
            end
            plot(W{ii}, 'LineWidth', 1, 'Color', c)
        %    plot(A{ii}, '*', 'LineWidth', 1, 'Color', c)
        end
    end

    % Unit circle boundary
    plot(xc, yc, 'Color', s.axiscolor, 'LineWidth', 1);

    srg_apply_style(s);
    xlim([-1 1])
    ylim([-1 1])

end

%--------------------------------------------------------------------------
function fill_bk_slice(W_slice, c, face_alpha)
%FILL_BK_SLICE  Fill a single W boundary in the BK disk.
    pts = W_slice(:);
    pts = pts(isfinite(pts));

    if length(pts) < 3
        return;
    end

    % Remove consecutive duplicates
    tol = max(1e-15, max(abs(pts)) * 1e-10);
    keep = [true; abs(diff(pts)) > tol];
    pts = pts(keep);

    if length(pts) < 3
        return;
    end

    c_pale = c + (1 - c)*0.45;

    ps = polyshape(real(pts), imag(pts), 'Simplify', false);
    if ps.NumRegions > 0
        plot(ps, 'FaceColor', c_pale, 'FaceAlpha', face_alpha, ...
             'EdgeColor', 'none');
    end
end

%--------------------------------------------------------------------------
function c = resolve_color(color, s)
%RESOLVE_COLOR  Always returns a plain 1x3 double RGB triplet.
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