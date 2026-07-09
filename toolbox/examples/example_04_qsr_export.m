%% example_04_qsr_export.m
%
% Example 4 — QSR multiplier extraction and publication export
% ============================================================
%
% PURPOSE
%   (1) srg_qsr_multiplier — solver-free (Q,S,R) multiplier from SRG slices.
%   (2) srg_bounding_circle — minimum enclosing circle of an SRG slice.
%   (3) srg_export — publication export.
%       Primary paths that always work:
%           .png   via getframe (no print driver needed)
%           .tikz  via matlab2tikz (vector, compile with pdflatex)
%       Best-effort (requires a working MATLAB print driver):
%           .pdf / .eps
%
% DEPENDENCIES
%   srg_compute, srg_qsr_multiplier, srg_bounding_circle,
%   srg_plot_gauss, srg_export
%   (matlab2tikz needed for the .tikz branch of Part 3)

clear; close all; clc;

rangemin  = -3;
rangemax  =  3;
estpoints = 200;

% Include \bg in the path
setup
%% =========================================================================
%  PART 1 — QSR multiplier of a MIMO 2x2 plant
% =========================================================================

G = [tf([1 2], [1 3 2])   tf([0.1], [1 1]);
     tf([0.05],[1 2])      tf([1],   [1 1])];

[~, ~, gplus, gminus, rawdata] = srg_compute(G, rangemin, rangemax, estpoints, 2000);

out = srg_qsr_multiplier(gplus, rawdata,'MaxOrder', 10, 'SKIter', 20, 'Plot', true);

fprintf('\nQSR multiplier summary\n');
fprintf('  S(s) order            : %d\n',   out.S_ord);
fprintf('  R(s) order            : %d\n',   out.R_ord);
fprintf('  R overbound factor k  : %.4g\n', out.R_scale);
fprintf('  pos-neg (r>|c|) holds : %d / %d frequencies\n',nnz(out.pos_neg), numel(out.pos_neg));

Pi = out.Pi; 
disp('Fitted S(s):');   disp(out.S_tf)
disp('Fitted R(s):');   disp(out.R_tf)




%% =========================================================================
%  PART 2 — Minimum enclosing circle of the widest SRG slice
% =========================================================================

G1=freqresp(G, 10^(-3), 'Hz');
figure
hold on;
[c, r] = srg_bounding_circle( gplus{1}(:), 'Plot', true);
srg_plot_gauss(G1,{gplus{1}},{gminus{1}},1,estpoints)
grid minor;
xlabel('Re'); ylabel('Im');
title(sprintf('SRG slice at %.3g Hz  (c=%.3f, r=%.3f)', 10^(-3), c, r));
legend({ 'min circle', 'centre','SRG'},'Location', 'best','Box','off','NumColumns',3);
xlim([0.85,1.15])
ylim([-0.15,0.15])

%% =========================================================================
%  PART 3 — 3-D surface SRG + per-frequency enclosing circle  (plant G)
% =========================================================================
fig3=figure;
srg_plot_3d_surface(rawdata,gplus, 0, 0,"FaceAlpha",0.9,"EdgeAlpha",1);   % surface SRG
hold on;

% Overlay minimum enclosing circle at every frequency
n=301;
th = linspace(0, 2*pi, n);
for i = 1:estpoints
    [c_i, r_i] = srg_bounding_circle(gplus{i}(:), 'Plot', false);
    xc = real(c_i) + r_i*cos(th);
    yc = imag(c_i) + r_i*sin(th);
    plot3(xc, yc, rawdata.freq(i)*ones(1,n), 'r-', 'LineWidth', 0.8);
end

set(gca, 'ZScale', 'log');
xlabel('Re'); ylabel('Im'); zlabel('Frequency (Hz)');


% =========================================================================
%  PART 3 — Publication export of a SISO SRG
% =========================================================================



% --- PNG (getframe-based: no print driver required) ----------------------
figure(fig3);
srg_export('example_04_srg.png', 'Width', 17, 'Height', 14.0);
fprintf('Wrote example_04_srg.png\n');

% --- TikZ (vector path to PDF: compile with pdflatex) -------------------
if exist('matlab2tikz', 'file')
    figure(fig3);
    srg_export('example_04_srg.tikz', 'Width', 8.5, 'Height', 7.0);
    fprintf('Wrote example_04_srg.tikz  (run: pdflatex example_04_srg.tikz)\n');
else
    fprintf('matlab2tikz not on path — skipping TikZ export.\n');
    fprintf('Install from https://github.com/matlab2tikz/matlab2tikz\n');
end

%% --- PDF (best-effort: requires a working MATLAB print driver) -----------
figure(fig3);
try
    srg_export('example_04_srg.pdf', 'Width', 17, 'Height', 14);
    fprintf('Wrote example_04_srg.pdf\n');
catch e
    fprintf('[info] PDF export unavailable on this system: %s\n', e.message);
    fprintf('       Use the PNG or TikZ output above instead.\n');
end