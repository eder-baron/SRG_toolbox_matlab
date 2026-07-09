function srg_export(filename, options)
%SRG_EXPORT  Export current figure in publication-ready format.
%
%   SRG_EXPORT(FILENAME) exports the current figure. The file extension
%   selects the backend:
%
%       .png   -- raster at specified DPI via getframe (no print driver needed)
%       .tikz  -- TikZ/PGFPlots via matlab2tikz (vector; compile with pdflatex)
%       .pdf   -- best-effort vector PDF (requires a working print driver)
%       .eps   -- best-effort vector EPS (requires a working print driver)
%
%   RECOMMENDED WORKFLOW for publication figures:
%       srg_export('fig.png',  'DPI', 300)      % always works
%       srg_export('fig.tikz')                   % vector via LaTeX (needs matlab2tikz)
%
%   Name-Value Arguments:
%       Fig        - Figure handle (default: gcf, evaluated at call time)
%       Width      - Width  in cm (default: 8.5 — one IEEE column)
%       Height     - Height in cm (default: 6.5)
%       DPI        - Raster resolution (default: 300, PNG only)
%       ApplyStyle - Re-apply publication theme before export (default: true)
%       FontSize   - Override font size; [] = theme default (default: [])
%
%   See also SRG_STYLE, SRG_APPLY_STYLE, SRG_PLOT_GAUSS

    arguments
        filename            (1,1) string
        options.Fig                     = []    % [] resolved to gcf below
        options.Width       (1,1) double = 8.5
        options.Height      (1,1) double = 6.5
        options.DPI         (1,1) double = 300
        options.ApplyStyle  (1,1) logical = true
        options.FontSize                = []
    end

    % ------------------------------------------------------------------
    %  Resolve figure: gcf is evaluated HERE (function body), not in the
    %  arguments block, to guarantee we get the runtime current figure.
    % ------------------------------------------------------------------
    if isempty(options.Fig)
        fig = gcf;
    else
        fig = options.Fig;
    end

    if ~isgraphics(fig, 'figure') || ~isvalid(fig)
        error('srg_export:invalidFig', ...
            ['No valid figure found.\n' ...
             'Call figure(fig) before srg_export, or pass ''Fig'', fig.']);
    end

    w = options.Width;
    h = options.Height;

    % ------------------------------------------------------------------
    %  Re-apply publication style to all axes and legends in the figure.
    % ------------------------------------------------------------------
    if options.ApplyStyle
        s = srg_style('Theme', 'publication');
        if ~isempty(options.FontSize)
            s.fontsize = options.FontSize;
        end
        for ax = findall(fig, 'Type', 'axes')'
            set(ax, 'FontSize', s.fontsize, 'FontName', s.fontname, ...
                    'TickLabelInterpreter', s.interpreter, 'Color', s.bgcolor);
            ax.XLabel.Interpreter = s.interpreter;
            ax.YLabel.Interpreter = s.interpreter;
            ax.Title.Interpreter  = s.interpreter;
            ax.XLabel.FontSize    = s.fontsize;
            ax.YLabel.FontSize    = s.fontsize;
        end
        for leg = findall(fig, 'Type', 'legend')'
            leg.Interpreter = s.interpreter;
            leg.FontSize    = s.fontsize;
            leg.FontName    = s.fontname;
        end
    end

    % ------------------------------------------------------------------
    %  Resize figure to the requested canvas (centimetres).
    % ------------------------------------------------------------------
    set(fig, 'Units', 'centimeters');
    set(fig, 'Position', [fig.Position(1:2), w, h]);

    % ------------------------------------------------------------------
    %  Dispatch by extension.
    % ------------------------------------------------------------------
    % Convert string -> char before fileparts so [fname,ext] concatenates
    % correctly; MATLAB string type makes [fname,ext] a 1x2 array, not a path.
    filename_c = char(filename);
    [fdir, fname, ext] = fileparts(filename_c);
    if isempty(fdir), fdir = '.'; end
    fullpath = fullfile(fdir, [fname, ext]);

    switch lower(ext)

        case '.png'
            srg_export_png(fig, fullpath, options.DPI);

        case {'.pdf', '.eps'}
            srg_export_vector(fig, fullpath, lower(ext));

        case '.tikz'
            srg_export_tikz(fig, fullpath, w, h);

        otherwise
            error('srg_export:unknownFormat', ...
                'Unsupported extension "%s". Use .png, .tikz, .pdf, or .eps.', ext);
    end

end


% ======================================================================
%  PNG  —  exportgraphics primary (clean crop); getframe fallback
% ======================================================================
function srg_export_png(fig, fullpath, target_dpi)
%SRGE_EXPORT_PNG  High-DPI raster export.
%   Tries exportgraphics first — it crops to the plot content and
%   handles DPI correctly without capturing figure chrome.
%   Falls back to a scaled getframe capture if exportgraphics fails.

    figure(fig);
    drawnow;

    % --- Primary: exportgraphics (same code path that works for PDF) ----
    try
        exportgraphics(gcf, fullpath, 'Resolution', target_dpi);
        fprintf('  srg_export: wrote %s  (~%d DPI)\n', fullpath, round(target_dpi));
        return;
    catch
    end

    % --- Fallback: getframe upscaling -----------------------------------
    %  getframe(fig) captures the full figure canvas including the gray
    %  background around the axes.  Scale the figure up first so the
    %  output is approximately target_dpi.
    screen_dpi = get(groot, 'ScreenPixelsPerInch');
    if isempty(screen_dpi) || screen_dpi <= 0
        screen_dpi = 96*5;
    end
    scale = target_dpi / screen_dpi;

    orig_units   = fig.Units;
    fig.Units    = 'pixels';
    orig_pos     = fig.Position;
    fig.Position = [orig_pos(1:2), round(orig_pos(3:4) * scale)];
    drawnow;

    frame = getframe(fig);

    fig.Position = orig_pos;
    fig.Units    = orig_units;
    drawnow;

    imwrite(frame.cdata, fullpath);
    fprintf('  srg_export: wrote %s  (%d x %d px, getframe fallback)\n', ...
        fullpath, size(frame.cdata, 2), size(frame.cdata, 1));
end


% ======================================================================
%  PDF / EPS  —  best-effort; clear message if unavailable
% ======================================================================
function srg_export_vector(fig, fullpath, ext)
%SRGE_EXPORT_VECTOR  Attempt vector export; fail clearly if print broken.

    is_pdf = strcmp(ext, '.pdf');
    fmt    = '-dpdf';
    if ~is_pdf, fmt = '-depsc'; end

    % Make sure the figure is the active current figure before any export call
    figure(fig);
    drawnow;

    % --- Attempt 1: exportgraphics (R2020a+, different stack from print) --
    try
        exportgraphics(gcf, fullpath, 'ContentType', 'vector');
        return;
    catch; end

    % --- Attempt 2: exportgraphics without ContentType -------------------
    try
        exportgraphics(gcf, fullpath);
        return;
    catch; end

    % --- Attempt 3: print without renderer flag --------------------------
    try
        print(gcf, fullpath, fmt);
        return;
    catch; end

    % --- Attempt 4: print with painters ----------------------------------
    try
        print(gcf, fullpath, fmt, '-painters');
        return;
    catch; end

    % All methods failed — give a clear message with working alternatives
    error('srg_export:vectorUnavailable', ...
        ['\n' ...
         'Vector export (%s) failed on this system.\n' ...
         'This usually means the MATLAB print driver is not configured.\n\n' ...
         'Working alternatives:\n' ...
         '  srg_export(''%s'', ''DPI'', 300)          %% high-res PNG\n' ...
         '  srg_export(''%s'', ''Width'', %.4g, ''Height'', %.4g)  %% TikZ -> pdflatex\n'], ...
        ext, ...
        strrep(fullpath, ext, '.png'), ...
        strrep(fullpath, ext, '.tikz'), ...
        fig.Position(3), fig.Position(4));
end


% ======================================================================
%  TikZ  —  matlab2tikz, ready for \input{} in LaTeX
% ======================================================================
function srg_export_tikz(fig, fullpath, w, h)
%SRGE_EXPORT_TIKZ  TikZ/PGFPlots output via matlab2tikz.

    if ~exist('matlab2tikz', 'file')
        error('srg_export:noMatlab2Tikz', ...
            ['matlab2tikz is not on the MATLAB path.\n' ...
             'Install from: https://github.com/matlab2tikz/matlab2tikz\n' ...
             'Then add to path: addpath(''/path/to/matlab2tikz/src'')']);
    end

    figure(fig);
    drawnow;

    matlab2tikz(fullpath, ...
        'figurehandle',    gcf, ...
        'width',           sprintf('%.4fcm', w), ...
        'height',          sprintf('%.4fcm', h), ...
        'strictFontSize',  false, ...
        'parseStrings',    true,  ...
        'showInfo',        false, ...
        'checkForUpdates', false);

    fprintf('  srg_export: wrote %s  (compile with pdflatex)\n', fullpath);
end