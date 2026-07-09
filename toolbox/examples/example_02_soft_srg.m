%% example_02_soft_srg.m
%
% Example 2 — Soft SRG (hyperbolic convex hull)
% ==============================================
%
% PURPOSE
%   Demonstrates the soft SRG, which is the hyperbolic convex hull of
%   the frequency-wise BK slices over all frequencies.  For stable,
%   minimum-phase systems it coincides with the operator-theoretic SRG.
%   Because Euclidean convex hull = hyperbolic convex hull in the BK
%   model, srg_scaled_graph operates directly on the W cell arrays.
%
% CASES
%   Loop A  –  SISO.    Two stable systems in negative feedback.
%   Loop B  –  MIMO 2x2.  Two stable systems in negative feedback.
%
% PIPELINE (per case)
%   1.  Define G1, G2 and confirm closed-loop stability (feedback / isstable)
%   2.  Compute frequency-wise SRG of each partner           (srg_compute)
%   3.  Compute soft SRG = hyperbolic convex hull            (srg_scaled_graph)
%   4.  Plot overlaid soft SRGs and compare with slices      (srg_plot_compare)
%
% DEPENDENCIES
%   srg_compute, srg_scaled_graph, srg_plot_gauss, srg_plot_beltrami,
%   srg_plot_compare

clear; close all; clc;
% Include \bg in the path
setup
% ── Shared sweep parameters ───────────────────────────────────────────────
rangemin  = -3;    % log10 of minimum frequency [Hz]
rangemax  =  3;    % log10 of maximum frequency [Hz]
estpoints = 100;   % number of frequency points

% points: angular resolution of the BK numerical range
%   SISO  → frequency slice is a scalar point, set points = 1
%   MIMO  → frequency slice is a curve, use points ~ 1000
points_siso = 1;
points_mimo = 2000;


% =========================================================================
%  LOOP A  —  SISO
% =========================================================================

G1_siso = tf([1 2], [2 3 2]);        % plant      — stable, zero at -2
G2_siso = tf([2 1], [2 1 3]);        % controller — stable

T_siso = feedback(G1_siso, G2_siso);
fprintf('--- Loop A (SISO):  G1 stable=%d, G2 stable=%d, T stable=%d\n', ...
        isstable(G1_siso), isstable(G2_siso), isstable(T_siso));

[W1, A1, gp1, gm1, ~] = srg_compute(G1_siso,       rangemin, rangemax, estpoints, points_siso);
[W2, A2, gp2, gm2, ~] = srg_compute(-inv(G2_siso),  rangemin, rangemax, estpoints, points_siso);

% Soft SRG = hyperbolic convex hull of all frequency slices
[sg1_plus, sg1_minus, h1] = srg_scaled_graph(W1);
[sg2_plus, sg2_minus, h2] = srg_scaled_graph(W2);

% -- 2D Gauss plane (soft SRG) ---------------------------------------------
figure('Name', 'Loop A — SISO, Gauss-plane soft SRGs');
srg_plot_gauss(G1_siso, sg1_plus, sg1_minus, 1, estpoints); hold on;
srg_plot_gauss(G2_siso, sg2_plus, sg2_minus, 2, estpoints);
axis equal; grid minor;
xlabel('Re'); ylabel('Im');
xlim([-4, 4]); ylim([-4, 4]);
title('Loop A: SISO — Gauss-plane soft SRGs');

% -- Overlaid soft SRGs ----------------------------------------------------
figure('Name', 'Loop A — SISO, soft SRGs overlaid');
srg_plot_compare(gp1, gp2, 'FaceAlpha', 0.7);
axis equal; grid minor;
xlabel('Re'); ylabel('Im');
xlim([-4, 4]); ylim([-4, 4]);
title('Loop A: SISO — soft SRG of plant vs controller');

% -- Freq-wise slices with soft hull overlaid ------------------------------
figure('Name', 'Loop A — SISO, slices vs hull');
srg_plot_gauss(G1_siso, sg1_plus, sg1_minus, 1, estpoints); hold on;
srg_plot_gauss(G2_siso, sg2_plus, sg2_minus, 2, estpoints);
srg_plot_compare(sg1_plus, sg2_plus, 'Fill', true, 'Color1', 1, 'Color2', 2);
axis equal; grid minor;
xlabel('Re'); ylabel('Im');
xlim([-4, 4]); ylim([-4, 4]);
title('Loop A: SISO — frequency slices with soft hull');


%% =========================================================================
%  LOOP B  —  MIMO 2x2
% =========================================================================

G1_mimo = [tf([1 2], [1 3 2])   tf([0.1],  [1 1]);
           tf([0.05], [1 2])     tf([1],    [1 1])];   % plant

G2_mimo = [tf([1 1], [1 5])     tf([1],    [3 2 2]);
           tf([1],   [1 2])    -tf([0.5],  [1 2])];    % controller

T_mimo = feedback(G1_mimo, G2_mimo);
fprintf('--- Loop B (MIMO 2x2): G1 stable=%d, G2 stable=%d, T stable=%d\n', ...
        isstable(G1_mimo), isstable(G2_mimo), isstable(T_mimo));

[W1, ~, gp1, gm1, ~] = srg_compute(G1_mimo,       rangemin, rangemax, estpoints, points_mimo);
[W2, ~, gp2, gm2, ~] = srg_compute(-inv(G2_mimo),  rangemin, rangemax, estpoints, points_mimo);

%% Soft SRG = hyperbolic convex hull of all frequency slices
[sg1_plus, ~, ~] = srg_scaled_graph(W1);
[sg2_plus, ~, ~] = srg_scaled_graph(W2);


%% -- Overlaid soft SRGs ----------------------------------------------------
figure('Name', 'Loop B — MIMO, soft SRGs overlaid');
srg_plot_compare(sg1_plus, sg2_plus, 'FaceAlpha', 0.8, 'Fill', true);
axis equal; grid minor;
xlabel('Re'); ylabel('Im');
xlim([-5, 5]); ylim([-5, 5]);
title('Loop B: MIMO 2x2 — soft SRG of plant vs controller');

%% -- Freq-wise slices vs soft hull (G1) ------------------------------------
figure('Name', 'Loop B — MIMO, G1 slices vs hull');
srg_plot_compare(gp1, sg1_plus, 'Names', ["freq-wise", "soft hull"]);
xlabel('Re'); ylabel('Im');
xlim([-.5, 1.5]); ylim([-1, 1]);
title('Loop B: MIMO 2x2 — G_1 frequency slices vs soft hull');

%% -- Freq-wise slices vs soft hull (G2) ------------------------------------
figure('Name', 'Loop B — MIMO, G2 slices vs hull');
srg_plot_compare(gp2, sg2_plus, 'Names', ["freq-wise", "soft hull"]);
axis equal; grid minor;
xlabel('Re'); ylabel('Im');
xlim([-5, 5]); ylim([-5, 5]);
title('Loop B: MIMO 2x2 — G_2 frequency slices vs soft hull');