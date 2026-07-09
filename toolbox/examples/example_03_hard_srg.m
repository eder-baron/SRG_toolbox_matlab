%% example_03_hard_srg.m
%
% Example 3 — Hard SRG for MIMO feedback pairs
% =============================================
%
% PURPOSE
%   Demonstrates srg_hard for MIMO LTI feedback interconnections where
%   one or both plants have unstable poles.  The soft (frequency-wise)
%   SRG cannot represent unbounded operators correctly; the hard SRG
%   handles this by capping the outer radius at a finite value (InfVal).
%
% CASES
%   Case A  –  MIMO 2x2.  G1 has an unstable pole.
%              Closed loop evaluated as feedback(G1, G2).
%   Case B  –  MIMO 3x3.  G2 has unstable poles and an integrator.
%              Closed loop evaluated as feedback(G2, G1).
%
% PIPELINE (per case)
%   1.  Define G1, G2 and confirm closed-loop stability (feedback / isstable)
%   2.  Set alpha/phi/frequency sweeps and InfVal cap
%   3.  Compute hard SRG of each partner                     (srg_hard)
%   4.  Plot overlaid filled regions                         (fill)
%
% DEPENDENCIES
%   srg_hard
 
clear; close all; clc;
 
 % Include \bg in the path
setup
% =========================================================================
%  CASE A  —  MIMO 2x2
% =========================================================================
 
s = tf('s');
 
% G1: 2x2 plant — entry (1,1) has an unstable pole at s = +1
G1 = [(s+7)/(s-1)         (s-5)/(s+2)^2; ...
       1/(s+4)^3           s/((s+3))^2  ] + 4*eye(2);
 
% G2: 2x2 controller — entry (2,1) has an unstable pole at s = +3
G2 = [(s+2)/((s+1)*(s+3))  (s+5)/(s+2)             ; ...
      (s+2)/(s-3)           (s+3)/((s+2)*(s+4))^2  ];
 
T = feedback(G1, G2);
fprintf('--- Case A (MIMO 2x2):\n');
fprintf('    G1 stable                   : %d\n', isstable(G1));
fprintf('    G2 stable                   : %d\n', isstable(G2));
fprintf('    Closed-loop feedback(G1,G2) : %d\n\n', isstable(T));
 
% ── Hard-SRG parameters ───────────────────────────────────────────────────
alphas = linspace(-10, 10, 100);   % alpha sweep (real-axis parameter)
phis   = linspace(0, pi,  100);    % phi sweep   (angle parameter)
freqs  = logspace(-3, 3,  500);     % frequency grid [rad/s]
infval = 50;                        % finite cap replacing Inf outer radius
limits = 6;                         % axis half-width for plotting [±limits]
 
fprintf('Computing hard SRG of G1 ...\n');
[srg1, ~] = srg_hard(G1,        alphas, phis, freqs, 'InfVal', infval);
 
fprintf('Computing hard SRG of -inv(G2) ...\n');
[srg2, ~] = srg_hard(-inv(G2),  alphas, phis, freqs, 'InfVal', infval);
 
% -- Plot ------------------------------------------------------------------
figure('Name', 'Case A — MIMO 2x2, hard SRGs');
fill(real(srg1), imag(srg1), [0.85 0.32 0.10], ...
     'FaceAlpha', 0.35, 'EdgeColor', 'none', 'DisplayName', 'G_1');
hold on;
fill(real(srg2), imag(srg2), [0.00 0.45 0.74], ...
     'FaceAlpha', 0.50, 'EdgeColor', 'none', 'DisplayName', '-inv(G_2)');
yline(0, 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
xline(0, 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
axis equal; grid minor;
xlabel('Re'); ylabel('Im');
xlim([-limits, limits]); ylim([-limits, limits]);
legend('Location', 'best');
title('Case A: MIMO 2x2 — hard SRGs of the feedback partners');
 
 
% =========================================================================
%  CASE B  —  MIMO 3x3
% =========================================================================
 
% G1: 3x3 stable plant
G1 = [(s+1)/(s+5)    (s+1)/(s+6)    (s+10)/(s+15); ...
      (s+5)/(s+11)   (s+10)/(s+15)  (2)/(s+1)    ; ...
      (5)/(s+1)      (1)/(s+6)      (3)/(s+1)    ] ...
     + [0.1  0   1 ; ...
        1    0.1 0 ; ...
        1.5  0  -1];
 
% G2: 3x3 plant — unstable pole at s = +1 in entry (1,1),
%                  integrator (1/s) in entry (3,3)
G2 = 15 * [(s+1)/(s-1)   (s-5)/(s+2)^2   (1)/(s+10); ...
            1/(s+4)^3    s/((s+3)^2)      (1)/(s+12); ...
            (2)/(s+1)    (1)/(s+14)        1/s      ] ...
     + 6*eye(3);
 
T = feedback(G2, G1);
fprintf('--- Case B (MIMO 3x3):\n');
fprintf('    G1 stable                   : %d\n', isstable(G1));
fprintf('    G2 stable                   : %d\n', isstable(G2));
fprintf('    Closed-loop feedback(G2,G1) : %d\n\n', isstable(T));
 
% ── Hard-SRG parameters ───────────────────────────────────────────────────
alphas = linspace(-10, 10, 1000);
phis   = linspace(0, pi,  1000);
freqs  = logspace(-3, 3,  200);     % coarser grid for speed (3x3 is heavier)
infval = 50;
limits = 10;
 
% G1 enters as the feedback partner, so its negated inverse is passed
fprintf('Computing hard SRG of -inv(G1) ...\n');
[srg1, ~] = srg_hard(-inv(G1), alphas, phis, freqs, 'InfVal', infval);
 
fprintf('Computing hard SRG of G2 ...\n');
[srg2, ~] = srg_hard(G2,       alphas, phis, freqs, 'InfVal', infval);
 
% -- Plot ------------------------------------------------------------------
figure('Name', 'Case B — MIMO 3x3, hard SRGs');
fill(real(srg1), imag(srg1), [0.85 0.32 0.10], ...
     'FaceAlpha', 0.35, 'EdgeColor', 'none', 'DisplayName', '-inv(G_1)');
hold on;
fill(real(srg2), imag(srg2), [0.00 0.45 0.74], ...
     'FaceAlpha', 0.50, 'EdgeColor', 'none', 'DisplayName', 'G_2');
yline(0, 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
xline(0, 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
axis equal; grid minor;
xlabel('Re'); ylabel('Im');
xlim([-limits, limits]); ylim([-limits, limits]);
legend('Location', 'best');
title('Case B: MIMO 3x3 — hard SRGs of the feedback partners');
 