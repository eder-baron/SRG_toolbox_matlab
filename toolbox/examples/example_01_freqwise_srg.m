%% example_01_freqwise_srg.m
%
% Example 1 — Frequency-wise SRG
% ===============================
%
% PURPOSE
%   Computes and plots the frequency-wise SRG of the feedback partners for
%   three interconnections of increasing complexity.  Each frequency slice
%   is a single BK-disk point (SISO) or a numerical-range arc (MIMO).
%
% CASES
%   Loop A  –  SISO.    Two stable systems in negative feedback.
%   Loop B  –  MIMO 2x2.  Two stable systems in negative feedback.
%   Loop C  –  MIMO 2x2.  Two stable systems; closed loop is unstable.
%
% PIPELINE (per case)
%   1.  Define G1, G2 and confirm closed-loop stability (feedback / isstable)
%   2.  Compute frequency-wise SRG of each partner           (srg_compute)
%   3.  Plot: 2D Gauss, Beltrami-Klein disk, 3D comparison   (srg_plot_*)
%   4.  Compute per-frequency stability margin                (srg_stability_margin)
%
% DEPENDENCIES
%   srg_compute, srg_plot_gauss, srg_plot_beltrami,
%   srg_plot_3d_compare, srg_plot_3d_beltrami_compare,
%   srg_stability_margin

clear; close all; clc;
% Include \bg in the path
setup
% ── Shared sweep parameters ───────────────────────────────────────────────
rangemin  = -3;    % log10 of minimum frequency [Hz]
rangemax  =  3;    % log10 of maximum frequency [Hz]
estpoints = 500;   % number of frequency points

% points: angular resolution of the BK numerical range
%   SISO  → frequency slice is a scalar point, set points = 1
%   MIMO  → frequency slice is a curve, use points ~ 500
points_siso = 1;
points_mimo = 500;


%% =========================================================================
%  LOOP A  —  SISO
% =========================================================================

G1_siso = tf([1 2], [1 3 2]);        % plant      — stable, zero at -2
G2_siso = tf([1 1], [1 5]);          % controller — stable

T_siso = feedback(G1_siso, G2_siso);
fprintf('--- Loop A (SISO):  G1 stable=%d, G2 stable=%d, T stable=%d\n', ...
        isstable(G1_siso), isstable(G2_siso), isstable(T_siso));

[W1, A1, gp1, gm1, rd1] = srg_compute(G1_siso,      rangemin, rangemax, estpoints, points_siso);
[W2, A2, gp2, gm2, rd2] = srg_compute(-inv(G2_siso), rangemin, rangemax, estpoints, points_siso);

% -- 2D Gauss plane --------------------------------------------------------
figure('Name', 'Loop A — SISO, Gauss-plane SRGs');
srg_plot_gauss(G1_siso, gp1, gm1, 4, estpoints); hold on;
srg_plot_gauss(G2_siso, gp2, gm2, 5, estpoints);
axis equal; grid minor;
xlabel('Re'); ylabel('Im');
title('Loop A: SISO — G_1 (plant) and G_2 (controller)');

% -- 2D Beltrami-Klein disk ------------------------------------------------
figure('Name', 'Loop A — SISO, Beltrami-Klein');
srg_plot_beltrami(G1_siso, W1, A1, 4, estpoints); hold on;
srg_plot_beltrami(G2_siso, W2, A2, 5, estpoints);
grid minor;
xlabel('Re'); ylabel('Im');
title('Loop A: SISO — Beltrami-Klein disk');

% -- 3D Gauss-plane comparison ---------------------------------------------
figure('Name', 'Loop A — SISO, 3D comparison');
srg_plot_3d_compare(rd1, gp1, rd2, gp2, 0, 0, ...
    'Color1', 4, 'Color2', 5, 'FaceAlpha', 0.4);
axis tight;
title('Loop A: SISO — 3D SRGs of G_1 and G_2');

% -- 3D Beltrami-Klein -----------------------------------------------------
figure('Name', 'Loop A — SISO, 3D Beltrami-Klein');
srg_plot_3d_beltrami_compare(rd1, W1, rd2, W2, 0, 0, ...
    'Color1', 4, 'Color2', 5, 'FaceAlpha', 0.8, ...
    'Lighting', true, 'IntersectColor', 1, 'IntersectAlpha', 1, ...
    'ShowDisk', true);
title('Loop A: SISO — 3D Beltrami-Klein disk');

% -- Per-frequency stability margin ----------------------------------------
[m_A, freq_A] = srg_stability_margin(gp1, gp2, rd1.freq);
figure('Name', 'Loop A — SISO, stability margin');
semilogx(freq_A, m_A, 'LineWidth', 1.5);
grid minor;
xlabel('Frequency (Hz)'); ylabel('SRG separation');
title(sprintf('Loop A: SISO stability margin  (min = %.3g)', min(m_A)));

% -- Step response ---------------------------------------------------------
figure('Name', 'Loop A — SISO, step response');
step(T_siso);
title('Loop A: SISO — step response');


%% =========================================================================
%  LOOP B  —  MIMO 2x2  (stable closed loop)
% =========================================================================

G1_mimo = [tf([1 2], [1 3 2])   tf([0.1],  [1 1]);
           tf([0.05], [1 2])     tf([1],    [1 1])];   % plant

G2_mimo = [tf([1 1], [1 5])     tf([1],    [3 2 2]);
           tf([1],   [1 2])    -tf([0.5],  [1 2])];    % controller

T_mimo = feedback(G1_mimo, G2_mimo, -1);
fprintf('--- Loop B (MIMO 2x2): G1 stable=%d, G2 stable=%d, T stable=%d\n', ...
        isstable(G1_mimo), isstable(G2_mimo), isstable(T_mimo));

[W1, A1, gp1, gm1, rd1] = srg_compute(-inv(G1_mimo), rangemin, rangemax, estpoints, points_mimo);
[W2, A2, gp2, gm2, rd2] = srg_compute(G2_mimo,       rangemin, rangemax, estpoints, points_mimo);

% -- 2D Gauss plane --------------------------------------------------------
figure('Name', 'Loop B — MIMO, Gauss-plane SRGs');
srg_plot_gauss(G1_mimo, gp1, gm1, 4, estpoints, 'Fill', true); hold on;
srg_plot_gauss(G2_mimo, gp2, gm2, 5, estpoints, 'Fill', true);
axis equal; grid minor;
xlabel('Re'); ylabel('Im');
xlim([-3, 3]); ylim([-3, 3]);
title('Loop B: MIMO 2x2 — G_1 (plant) and G_2 (controller)');

% -- 2D Beltrami-Klein disk ------------------------------------------------
figure('Name', 'Loop B — MIMO, Beltrami-Klein');
srg_plot_beltrami(G1_mimo, W1, A1, 4, estpoints, 'Fill', true); hold on;
srg_plot_beltrami(G2_mimo, W2, A2, 5, estpoints, 'Fill', true);
grid minor;
xlabel('Re'); ylabel('Im');
title('Loop B: MIMO 2x2 — Beltrami-Klein disk');

% -- 3D Gauss-plane comparison ---------------------------------------------
figure('Name', 'Loop B — MIMO, 3D comparison');
srg_plot_3d_compare(rd1, gp1, rd2, gp2, 0, 0, ...
    'Color1', 4, 'Color2', 5, 'FaceAlpha', 0.8, 'EdgeAlpha', 0.6, ...
    'Lighting', false);
view(50, 25); axis tight;
xlim([-3, 3]); ylim([-3, 3]);
title('Loop B: MIMO 2x2 — 3D SRGs of G_1 and G_2');
%%
figure
srg_plot_3d_compare(rd1, gp1, rd2, gp2, 0, 0, ...
    'Color1', 4, 'Color2', 5, 'FaceAlpha', 0.8, 'EdgeAlpha', 0.6, ...
    'Lighting', false);
[sm,fr,x1,x2] = srg_stability_margin(gp1,gp2,rd1.freq);
srg_plot_witness_3d(x1,x2,fr,sm,'Stride',10,'LineWidth',2,'MarkerSize',0.1,'Color',1);
xlim([-3, 3]); ylim([-3, 3]);


%% -- 3D Beltrami-Klein -----------------------------------------------------
figure('Name', 'Loop B — MIMO, 3D Beltrami-Klein');
srg_plot_3d_beltrami_compare(rd1, W1, rd2, W2, 0, 0, ...
    'Color1', 4, 'Color2', 5, 'FaceAlpha', 0.8, ...
    'Lighting', true, 'IntersectColor', 1, 'IntersectAlpha', 1);
title('Loop B: MIMO 2x2 — 3D Beltrami-Klein disk');

% -- Per-frequency stability margin ----------------------------------------
[m_B, freq_B] = srg_stability_margin(gp1, gp2, rd1.freq);
figure('Name', 'Loop B — MIMO, stability margin');
semilogx(freq_B, m_B, 'LineWidth', 1.5);
grid minor;
xlabel('Frequency (Hz)'); ylabel('SRG separation');
ylim([0, 4]);
title(sprintf('Loop B: MIMO 2x2 stability margin  (min = %.3g)', min(m_B)));

% -- Step response ---------------------------------------------------------
figure('Name', 'Loop B — MIMO, step response');
step(T_mimo);
title('Loop B: MIMO 2x2 — step response');

% Homotopy
hd = srg_homotopy(G2_mimo, [1.0, 0.6, 0.3], -3, 3, 100, 64);
figure;
for k = 1:numel(hd)
srg_plot_3d_surface(hd(k).rawdata, hd(k).gplus, 0, 0,'Color', k,'FaceAlpha', 0.8, 'EdgeAlpha', 0.6, ...
    'Lighting', false);
hold on
end
srg_plot_3d_surface(rd1, gp1, 0, 0,'Color', 4,'FaceAlpha', 0.8, 'EdgeAlpha', 0.6, ...
    'Lighting', false);
xlim([-3, 3]); ylim([-3, 3]);

%% =========================================================================
%  LOOP C  —  MIMO 2x2  (unstable closed loop)
% =========================================================================

H1_mimo = [tf([1],    [1 1])   tf([1],   [1 2]);    % DC: [1.0  0.5]
           tf([1.5],  [1 3])   tf([4],   [1 4])];   % DC: [0.5  1.0]

H2_mimo = [tf([0.25], [1 5])   tf([18],  [1 6]);    % DC: [0.05  3.0]
           tf([21],   [1 7])   tf([0.4], [1 8])];   % DC: [3.0   0.05]

L_mimo = feedback(H1_mimo, H2_mimo, -1);
fprintf('--- Loop C (MIMO 2x2): G1 stable=%d, G2 stable=%d, T stable=%d\n', ...
        isstable(H1_mimo), isstable(H2_mimo), isstable(L_mimo));

[W1, A1, gp1, gm1, rd1] = srg_compute(H1_mimo,       rangemin, rangemax, estpoints, points_mimo);
[W2, A2, gp2, gm2, rd2] = srg_compute(-inv(H2_mimo),  rangemin, rangemax, estpoints, points_mimo);

% -- 2D Gauss plane --------------------------------------------------------
figure('Name', 'Loop C — MIMO, Gauss-plane SRGs');
srg_plot_gauss(H1_mimo, gp1, gm1, 4, estpoints, 'Fill', true); hold on;
srg_plot_gauss(H2_mimo, gp2, gm2, 5, estpoints, 'Fill', true);
axis equal; grid minor;
xlabel('Re'); ylabel('Im');
xlim([-1, 1]); ylim([-1, 1]);
title('Loop C: MIMO 2x2 (unstable) — G_1 and G_2');

% -- 2D Beltrami-Klein disk ------------------------------------------------
figure('Name', 'Loop C — MIMO, Beltrami-Klein');
srg_plot_beltrami(H1_mimo, W1, A1, 4, estpoints, 'Fill', true); hold on;
srg_plot_beltrami(H2_mimo, W2, A2, 5, estpoints, 'Fill', true);
grid minor;
xlabel('Re'); ylabel('Im');
title('Loop C: MIMO 2x2 (unstable) — Beltrami-Klein disk');

% -- 3D Gauss-plane comparison ---------------------------------------------
figure('Name', 'Loop C — MIMO, 3D comparison');
srg_plot_3d_compare(rd1, gp1, rd2, gp2, 0, 0, ...
    'Color1', 4, 'Color2', 5, 'EdgeAlpha', 0.3, 'FaceAlpha', 0.8, ...
    'Lighting', true, 'IntersectColor', 1, 'IntersectAlpha', 1);
view(50, 25); axis tight;
xlim([-1, 1]); ylim([-1, 1]);
title('Loop C: MIMO 2x2 (unstable) — 3D SRGs of G_1 and G_2');

% -- 3D Beltrami-Klein -----------------------------------------------------
figure('Name', 'Loop C — MIMO, 3D Beltrami-Klein');
srg_plot_3d_beltrami_compare(rd1, W1, rd2, W2, 0, 0, ...
    'Color1', 4, 'Color2', 5, 'FaceAlpha', 0.8, ...
    'Lighting', true, 'IntersectColor', 1, 'IntersectAlpha', 1);
title('Loop C: MIMO 2x2 (unstable) — 3D Beltrami-Klein disk');

% -- Per-frequency stability margin ----------------------------------------
[m_C, freq_C] = srg_stability_margin(gp1, gp2, rd1.freq);
figure('Name', 'Loop C — MIMO, stability margin');
semilogx(freq_C, m_C, 'LineWidth', 1.5);
grid minor;
xlabel('Frequency (Hz)'); ylabel('SRG separation');
ylim([0, 3]);
title(sprintf('Loop C: MIMO 2x2 stability margin  (min = %.3g)', min(m_C)));
