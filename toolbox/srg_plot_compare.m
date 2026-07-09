function srg_plot_compare(varargin)
%SRG_PLOT_COMPARE  Overlay multiple SRGs with automatic styling.
%
%   SRG_PLOT_COMPARE(SRG1, SRG2, ...) plots multiple SRG boundaries as
%   filled regions with automatic color cycling, clean polyshape rendering,
%   and a legend.
%
%   Name-Value Arguments (must come after all SRG vectors):
%       Names     - String array of legend names  (default: "SRG 1", ...)
%       Theme     - Style theme                   (default: "default")
%       Mode      - "filled" (default) or "boundary"
%       FaceAlpha - Fill transparency override    (default: from theme)
%       Title     - Plot title                    (default: '')
%       Fill      - Fill interior of each SRG cell (default: false)
%                   When true and SRG inputs are cell arrays of cell arrays
%                   (i.e. gplus from srg_compute), fills each frequency
%                   slice.  For SISO inputs (each cell is a single point)
%                   fills the region swept by the full trajectory.
%                   When SRG inputs are plain complex vectors (e.g. from
%                   srg_soft), uses polyshape fill as before.
%
%   SISO detection: if every cell in a cell-array SRG contains a single
%   unique complex point the input is treated as SISO and the SRG is drawn
%   as a connected frequency-sweep trajectory rather than per-slice curves.
%
%   See also SRG_HARD, SRG_SCALED_GRAPH, SRG_PLOT_GAUSS, SRG_COMPUTE

    % ------------------------------------------------------------------
    % Split positional SRG arguments from name-value pairs
    % ------------------------------------------------------------------
    n_srg = 0;
    for kk = 1:nargin
        if isnumeric(varargin{kk}) || iscell(varargin{kk})
            n_srg = n_srg + 1;
        else
            break;
        end
    end

    srgs    = varargin(1:n_srg);
    nv_args = varargin(n_srg+1:end);

    % ------------------------------------------------------------------
    % Defaults
    % ------------------------------------------------------------------
    names = strings(1, n_srg);
    for kk = 1:n_srg
        names(kk) = sprintf("SRG %d", kk);
    end
    theme      = "default";
    mode       = "filled";
    face_alpha = [];
    plot_title = '';
    do_fill    = false;

    % ------------------------------------------------------------------
    % Parse name-value pairs
    % ------------------------------------------------------------------
    ii = 1;
    while ii <= length(nv_args)
        key = string(nv_args{ii});
        switch lower(key)
            case "names"
                names = string(nv_args{ii+1});  ii = ii + 2;
            case "theme"
                theme = string(nv_args{ii+1});  ii = ii + 2;
            case "mode"
                mode  = string(nv_args{ii+1});  ii = ii + 2;
            case "facealpha"
                face_alpha = nv_args{ii+1};     ii = ii + 2;
            case "title"
                plot_title = nv_args{ii+1};     ii = ii + 2;
            case "fill"
                do_fill = logical(nv_args{ii+1}); ii = ii + 2;
            otherwise
                ii = ii + 1;
        end
    end

    s = srg_style('Theme', theme);
    if isempty(face_alpha)
        face_alpha = s.facealpha;
    end

    % ------------------------------------------------------------------
    % Plot each SRG
    % ------------------------------------------------------------------
    hold on
    legend_handles = gobjects(n_srg, 1);

    for kk = 1:n_srg
        c      = s.color(kk);
        c_pale = c + (1 - c) * 0.45;
        data   = srgs{kk};

        is_cell_of_curves = iscell(data);

        % ==============================================================
        % Cell-array input (from srg_compute: gplus / W)
        % ==============================================================
        if is_cell_of_curves

            if is_siso_cell(data)
                % ----------------------------------------------------------
                % SISO branch — each cell is one unique point.
                % Collect trajectory across frequencies and plot as a curve.
                % ----------------------------------------------------------
                traj = collect_siso_trajectory(data);   % complex vector

                if do_fill
                    % Fill the region enclosed by upper (gplus) and lower
                    % (gminus = conj) trajectories as a single polyshape.
                    boundary_full = [traj; conj(flipud(traj))];
                    ps = polyshape(real(boundary_full), imag(boundary_full), ...
                                   'Simplify', false);
                    if ps.NumRegions > 0
                        h = plot(ps, 'FaceColor', c_pale, ...
                                 'FaceAlpha', face_alpha, ...
                                 'EdgeColor', c, 'LineWidth', s.linewidth);
                    else
                        % Fallback: just draw the trajectory lines
                        plot(real(traj), imag(traj), '-', ...
                             'Color', c, 'LineWidth', s.linewidth, ...
                             'HandleVisibility', 'off');
                        plot(real(traj), -imag(traj), '-', ...
                             'Color', c, 'LineWidth', s.linewidth, ...
                             'HandleVisibility', 'off');
                        h = fill(NaN, NaN, c_pale, 'FaceAlpha', face_alpha, ...
                                 'EdgeColor', c, 'LineWidth', s.linewidth);
                    end
                else
                    % Boundary only: upper + lower (conjugate) trajectories
                    plot(real(traj),  imag(traj), '-', ...
                         'Color', c, 'LineWidth', s.linewidth, ...
                         'HandleVisibility', 'off');
                    plot(real(traj), -imag(traj), '-', ...
                         'Color', c, 'LineWidth', s.linewidth, ...
                         'HandleVisibility', 'off');
                    h = plot(NaN, NaN, '-', 'Color', c, ...
                             'LineWidth', s.linewidth);
                end
                legend_handles(kk) = h;

            elseif do_fill
                % ----------------------------------------------------------
                % MIMO + fill — fill each frequency slice individually
                % ----------------------------------------------------------
                for jj = 1:length(data)
                    curve = clean_curve(data{jj});
                    if length(curve) < 3, continue; end

                    boundary_full = [curve; conj(flipud(curve))];
                    ps = polyshape(real(boundary_full), imag(boundary_full), ...
                                   'Simplify', false);
                    if ps.NumRegions > 0
                        plot(ps, 'FaceColor', c_pale, ...
                             'FaceAlpha', face_alpha, 'EdgeColor', 'none');
                    end
                end
                h = fill(NaN, NaN, c_pale, 'FaceAlpha', face_alpha, ...
                         'EdgeColor', c, 'LineWidth', s.linewidth);
                legend_handles(kk) = h;

            else
                % ----------------------------------------------------------
                % MIMO + boundary only — plot each slice as a curve
                % ----------------------------------------------------------
                for jj = 1:length(data)
                    curve = data{jj};
                    if isempty(curve), continue; end
                    plot(real(curve),  imag(curve), '-', ...
                         'Color', c, 'LineWidth', s.linewidth, ...
                         'HandleVisibility', 'off');
                    plot(real(curve), -imag(curve), '-', ...
                         'Color', c, 'LineWidth', s.linewidth, ...
                         'HandleVisibility', 'off');
                end
                h = plot(NaN, NaN, '-', 'Color', c, 'LineWidth', s.linewidth);
                legend_handles(kk) = h;
            end

        % ==============================================================
        % Plain complex vector input (from srg_soft / srg_hard)
        % ==============================================================
        elseif mode == "filled"
            ps = srg_to_polyshape(data);
            if ps.NumRegions > 0
                h = plot(ps, 'FaceColor', c, 'FaceAlpha', face_alpha, ...
                         'EdgeColor', c, 'EdgeAlpha', 0.6, 'LineWidth', 0.8);
            else
                [~, xs, ys] = srg_to_polyshape(data);
                h = fill(xs, ys, c, 'FaceAlpha', face_alpha, ...
                         'EdgeColor', c, 'EdgeAlpha', 0.6, 'LineWidth', 0.8);
            end
            legend_handles(kk) = h;

        else
            % boundary only
            [~, xs, ys] = srg_to_polyshape(data);
            xs = [xs; xs(1)];
            ys = [ys; ys(1)];
            h  = plot(xs, ys, '-', 'Color', c, 'LineWidth', s.linewidth);
            legend_handles(kk) = h;
        end
    end

    % ------------------------------------------------------------------
    % Apply style, legend, title
    % ------------------------------------------------------------------
    srg_apply_style(s);

    legend(legend_handles, names, 'Location', 'best', ...
           'FontSize', s.fontsize - 2, 'FontName', s.fontname, ...
           'Interpreter', s.interpreter, 'Box', 'off');

    if ~isempty(plot_title)
        title(plot_title, 'FontSize', s.fontsize + 2, ...
              'FontName', s.fontname, 'Interpreter', s.interpreter);
    end

end


% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function tf = is_siso_cell(gp)
%IS_SISO_CELL  True when every non-empty cell contains one unique point.
%   This is the signature of a SISO system processed through srg_compute:
%   gplus{i} holds one complex value (possibly repeated) per frequency.
    tol = 1e-8;
    tf  = true;
    for ii = 1:numel(gp)
        c = gp{ii};
        c = c(isfinite(c));
        if isempty(c), continue; end
        if max(abs(c - c(1))) > tol * max(abs(c(1)), 1)
            tf = false;
            return;
        end
    end
end


function traj = collect_siso_trajectory(gp)
%COLLECT_SISO_TRAJECTORY  Extract one representative point per cell.
    n    = numel(gp);
    traj = zeros(n, 1);
    for ii = 1:n
        c = gp{ii};
        c = c(isfinite(c));
        if isempty(c)
            traj(ii) = NaN;
        else
            traj(ii) = c(1);
        end
    end
    traj = traj(isfinite(traj));
end


function curve = clean_curve(raw)
%CLEAN_CURVE  Remove non-finite values and consecutive near-duplicates.
    curve = raw(isfinite(raw));
    curve = curve(:);
    if length(curve) < 2, return; end
    tol  = max(1e-15, max(abs(curve)) * 1e-10);
    keep = [true; abs(diff(curve)) > tol];
    curve = curve(keep);
end