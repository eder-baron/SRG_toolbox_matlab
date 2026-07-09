function srg_apply_style(s, options)
%SRG_APPLY_STYLE  Apply SRG style configuration to the current axes.
%
%   SRG_APPLY_STYLE(S) applies fonts, grid, background, and axis
%   crosshairs from the style struct S (from srg_style).
%
%   Name-Value Arguments:
%       Crosshairs - Draw x=0 and y=0 lines (default: true)
%       AxisEqual  - Set axis equal (default: true)
%       Labels     - Cell of {xlabel, ylabel} (default: {'Re','Im'})

    arguments
        s           struct
        options.Crosshairs (1,1) logical = true
        options.AxisEqual  (1,1) logical = true
        options.Labels     cell = {'Re', 'Im'}
    end

    ax = gca;

    % Font
    set(ax, 'FontSize', s.fontsize, 'FontName', s.fontname);

    % Background
    set(ax, 'Color', s.bgcolor);

    % Grid
    grid(ax, s.grid);

    % Labels
    if ~isempty(options.Labels)
        xlabel(options.Labels{1}, 'FontSize', s.fontsize, ...
               'FontName', s.fontname, 'Interpreter', s.interpreter);
        if length(options.Labels) >= 2
            ylabel(options.Labels{2}, 'FontSize', s.fontsize, ...
                   'FontName', s.fontname, 'Interpreter', s.interpreter);
        end
        if length(options.Labels) >= 3
            zlabel(options.Labels{3}, 'FontSize', s.fontsize, ...
                   'FontName', s.fontname, 'Interpreter', s.interpreter);
        end
    end

    % Axis equal
    if options.AxisEqual
        axis equal
    end

    % Crosshairs
    if options.Crosshairs
        hold on
        xline(0, 'Color', s.axiscolor, 'LineWidth', 0.5);
        yline(0, 'Color', s.axiscolor, 'LineWidth', 0.5);
    end

    % Box
    box on

end
