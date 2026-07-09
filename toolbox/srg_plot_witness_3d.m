function srg_plot_witness_3d(x1, x2, freq, stab_margin, options)
%SRG_PLOT_WITNESS_3D  Overlay minimum-distance witness segments on a 3D SRG plot.
%
%   SRG_PLOT_WITNESS_3D(X1, X2, FREQ, STAB_MARGIN) draws the single critical
%   witness segment (at the frequency of minimum stability margin) and marks
%   both endpoints.  Call after SRG_PLOT_3D_COMPARE on the same axes.
%
%   Inputs:
%       x1          - Complex vector, closest point on SRG 1 at each frequency
%       x2          - Complex vector, closest point on SRG 2 at each frequency
%       freq        - Frequency vector in Hz
%       stab_margin - Minimum distance at each frequency
%
%   Name-Value Arguments:
%       Theme      - Style theme (default: "default")
%       Mirror     - Also draw conjugate lower half-plane segments (default: true)
%       Color      - Segment color, palette index or RGB (default: 3)
%       LineWidth  - Segment line width (default: 2)
%       MarkerSize - Endpoint marker size (default: 10)
%       Stride     - Also draw every Stride-th background segment in thin grey.
%                    0 = critical only (default: 0)
%
%   Example:
%       [~,~,gp1,~,rd1] = srg_compute(G1,-2,3,200,32);
%       [~,~,gp2,~,rd2] = srg_compute(G2,-2,3,200,32);
%       figure;
%       srg_plot_3d_compare(rd1,gp1,rd2,gp2,0,0);
%       [sm,fr,x1,x2] = srg_stability_margin(gp1,gp2,rd1.freq);
%       srg_plot_witness_3d(x1,x2,fr,sm);
%       % or, to also show every 10th background segment:
%       srg_plot_witness_3d(x1,x2,fr,sm,'Stride',10);
%
%   See also SRG_STABILITY_MARGIN, SRG_PLOT_3D_COMPARE

    arguments
        x1            (:,1) double
        x2            (:,1) double
        freq          (:,1) double
        stab_margin   (:,1) double
        options.Theme      string  = "default"
        options.Mirror     (1,1) logical = true
        options.Color               = 8    % cb_orange: contrasts with blue(1) and purple(5)
        options.LineWidth  (1,1) double = 2
        options.MarkerSize (1,1) double = 10
        options.Stride     (1,1) double {mustBeNonnegative, mustBeInteger} = 0
    end

    s  = srg_style('Theme', options.Theme);
    c  = resolve_color(options.Color, s);
    lw = options.LineWidth;

    hold on;

    % -----------------------------------------------------------------------
    %  Optional background: every Stride-th segment in thin grey
    % -----------------------------------------------------------------------
    if options.Stride > 0
        idx_bg = 1 : options.Stride : numel(freq);
        bg = [0.7 0.7 0.7];

        for ii = idx_bg
            plot3([real(x1(ii)) real(x2(ii))], ...
                  [imag(x1(ii)) imag(x2(ii))], ...
                  [freq(ii)     freq(ii)], ...
                  '-', 'Color', bg, 'LineWidth', 0.5);
            if options.Mirror
                plot3([real(x1(ii)) real(x2(ii))], ...
                      [-imag(x1(ii)) -imag(x2(ii))], ...
                      [freq(ii)      freq(ii)], ...
                      '-', 'Color', bg, 'LineWidth', 0.5);
            end
        end
    end

    % -----------------------------------------------------------------------
    %  Critical segment: shortest witness segment only
    % -----------------------------------------------------------------------
    [sm_crit, i_crit] = min(stab_margin);
    fc  = freq(i_crit);
    p1c = x1(i_crit);
    p2c = x2(i_crit);

    plot3([real(p1c) real(p2c)], [imag(p1c) imag(p2c)], [fc fc], ...
          '-', 'Color', c, 'LineWidth', lw);
    plot3(real(p1c), imag(p1c), fc, 'o', ...
          'Color', c, 'MarkerFaceColor', c, 'MarkerSize', options.MarkerSize);
    plot3(real(p2c), imag(p2c), fc, 'o', ...
          'Color', c, 'MarkerFaceColor', c, 'MarkerSize', options.MarkerSize);

    if options.Mirror
        plot3([real(p1c) real(p2c)], [-imag(p1c) -imag(p2c)], [fc fc], ...
              '-', 'Color', c, 'LineWidth', lw);
        plot3(real(p1c), -imag(p1c), fc, 'o', ...
              'Color', c, 'MarkerFaceColor', c, 'MarkerSize', options.MarkerSize);
        plot3(real(p2c), -imag(p2c), fc, 'o', ...
              'Color', c, 'MarkerFaceColor', c, 'MarkerSize', options.MarkerSize);
    end

    fprintf('Minimum stability margin: %.4g  at  %.4g Hz\n', sm_crit, fc);

end

% --------------------------------------------------------------------------
function c = resolve_color(color, s)
%RESOLVE_COLOR  Palette index (via s.color function handle) or RGB passthrough.
    if isnumeric(color) && numel(color) == 3
        c = color(:)';
    else
        c = s.color(round(color));   % s.color is @(idx) s.palette(mod...) 
    end
end