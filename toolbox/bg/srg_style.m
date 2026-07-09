function s = srg_style(options)
%SRG_STYLE  Unified style configuration for SRG plots.
%
%   S = SRG_STYLE() returns a struct with default plot styling.
%   S = SRG_STYLE('Theme', 'light') or 'dark' or 'publication'.
%
%   The struct S contains:
%       s.palette     - Nx3 color matrix (curated for SRG overlays)
%       s.linewidth   - Default line width
%       s.fontsize    - Font size for labels
%       s.fontname    - Font name
%       s.facealpha   - Default fill transparency
%       s.interpreter - Text interpreter ('latex' or 'tex')
%       s.grid        - Grid style ('minor' or 'on')
%       s.axiscolor   - Color for axis crosshairs
%       s.bgcolor     - Axes background color
%
%   Usage in other functions:
%       s = srg_style();
%       c = s.palette(mod(idx-1, size(s.palette,1))+1, :);

    arguments
        options.Theme string {mustBeMember(options.Theme, ...
            ["default","light","dark","publication"])} = "default"
    end

    % CB Colorblind Safe palette (from IPE stylesheets)
    % Primary colors: good contrast, safe for all color vision types
    cb_blue   = [0.2667 0.4667 0.6667];
    cb_cyan   = [0.4000 0.8000 0.9333];
    cb_green  = [0.1333 0.5333 0.2000];
    cb_yellow = [0.8000 0.7333 0.2667];
    cb_red    = [0.9333 0.4000 0.4667];
    cb_purple = [0.6667 0.2000 0.4667];
    cb_grey   = [0.7333 0.7333 0.7333];
    cb_orange = [0.8500 0.3250 0.0980];

    % Pale variants: fills, backgrounds, secondary elements
    cb_pale_blue   = [0.7333 0.8000 0.9333];
    cb_pale_cyan   = [0.8000 0.9333 1.0000];
    cb_pale_green  = [0.8000 0.8667 0.6667];
    cb_pale_yellow = [0.9333 0.9333 0.7333];
    cb_pale_red    = [1.0000 0.8000 0.8000];
    cb_pale_grey   = [0.8667 0.8667 0.8667];

    % Dark variants: dark theme line colors, emphasis
    cb_dark_blue   = [0.1333 0.1333 0.3333];
    cb_dark_cyan   = [0.1333 0.3333 0.3333];
    cb_dark_green  = [0.1333 0.3333 0.1333];
    cb_dark_yellow = [0.4000 0.4000 0.2000];
    cb_dark_red    = [0.4000 0.2000 0.2000];
    cb_dark_grey   = [0.3333 0.3333 0.3333];

    % Light/default/publication: primary CB colors
    light_palette = [
        cb_blue
        cb_red
        cb_green
        cb_yellow
        cb_purple
        cb_cyan
        cb_grey
        cb_orange
    ];

    % Dark theme: use the pale variants (high visibility on dark bg)
    dark_palette = [
        cb_pale_blue
        cb_pale_red
        cb_cyan
        cb_pale_yellow
        cb_pale_green
        cb_pale_grey
        cb_grey
    ];

    % Store all CB colors for direct access
    s.cb.blue        = cb_blue;
    s.cb.cyan        = cb_cyan;
    s.cb.green       = cb_green;
    s.cb.yellow      = cb_yellow;
    s.cb.red         = cb_red;
    s.cb.purple      = cb_purple;
    s.cb.grey        = cb_grey;
    s.cb.orange      = cb_orange;
    s.cb.pale_blue   = cb_pale_blue;
    s.cb.pale_cyan   = cb_pale_cyan;
    s.cb.pale_green  = cb_pale_green;
    s.cb.pale_yellow = cb_pale_yellow;
    s.cb.pale_red    = cb_pale_red;
    s.cb.pale_grey   = cb_pale_grey;
    s.cb.dark_blue   = cb_dark_blue;
    s.cb.dark_cyan   = cb_dark_cyan;
    s.cb.dark_green  = cb_dark_green;
    s.cb.dark_yellow = cb_dark_yellow;
    s.cb.dark_red    = cb_dark_red;
    s.cb.dark_grey   = cb_dark_grey;

    switch options.Theme
        case "publication"
            s.palette     = light_palette;
            s.linewidth   = 1.5;
            s.fontsize    = 14;
            s.fontname    = 'Times New Roman';
            s.facealpha   = 0.35;
            s.interpreter = 'latex';
            s.grid        = 'on';
            s.axiscolor   = cb_dark_grey;
            s.bgcolor     = [1 1 1];

        case "dark"
            s.palette     = dark_palette;
            s.linewidth   = 1.8;
            s.fontsize    = 12;
            s.fontname    = 'Helvetica';
            s.facealpha   = 0.5;
            s.interpreter = 'tex';
            s.grid        = 'on';
            s.axiscolor   = cb_grey;
            s.bgcolor     = [0.15 0.15 0.18];

        case "light"
            s.palette     = light_palette;
            s.linewidth   = 1.5;
            s.fontsize    = 12;
            s.fontname    = 'Helvetica';
            s.facealpha   = 0.4;
            s.interpreter = 'tex';
            s.grid        = 'minor';
            s.axiscolor   = cb_dark_grey;
            s.bgcolor     = [1 1 1];

        otherwise % "default"
            s.palette     = light_palette;
            s.linewidth   = 1.5;
            s.fontsize    = 12;
            s.fontname    = 'Helvetica';
            s.facealpha   = 0.4;
            s.interpreter = 'tex';
            s.grid        = 'minor';
            s.axiscolor   = cb_dark_grey;
            s.bgcolor     = [1 1 1];
    end

    % Convenience: function to cycle palette
    % If this conflicts with the Symbolic Toolbox, use srg_get_color(s, N) instead.
    s.color = @(idx) s.palette(mod(idx-1, size(s.palette,1)) + 1, :);
end