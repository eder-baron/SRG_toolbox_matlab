function homotopy_data = srg_homotopy(G, tau_values, rangemin, rangemax, ...
                                       estpoints, points, options)
%SRG_HOMOTOPY  Frequency-wise SRG of τ·G for multiple τ values.
%
%   HOMOTOPY_DATA = SRG_HOMOTOPY(G, TAU_VALUES, RANGEMIN, RANGEMAX,
%   ESTPOINTS, POINTS) computes the full frequency-wise SRG of the
%   scaled system G_τ = τ·G at each τ in TAU_VALUES, by calling
%   SRG_COMPUTE once per τ.
%
%   Each entry of HOMOTOPY_DATA contains a complete freq-wise SRG
%   dataset for one τ — ready to feed individually into SRG_PLOT_3D,
%   SRG_PLOT_3D_SURFACE, SRG_PLOT_GAUSS, or SRG_PLOT_COMPARE.
%
%   Inputs:
%       G          - Transfer function (tf/ss object)
%       tau_values - Vector of τ values, typically in [0, 1]
%                    (e.g. [1.0, 0.6, 0.3] for three overlaid surfaces)
%       rangemin   - Minimum frequency exponent (log) or value (linear)
%       rangemax   - Maximum frequency exponent (log) or value (linear)
%       estpoints  - Number of frequency evaluation points
%       points     - Angular resolution for field of values
%
%   Name-Value Arguments:
%       FreqScale  - Frequency spacing: "log" (default) or "linear"
%
%   Output:
%       homotopy_data - Struct array of length numel(tau_values).
%                       Each element has fields:
%           .tau     - scalar τ value
%           .W       - cell array, BK numerical range boundaries
%           .A       - cell array, eigenvalues in BK plane
%           .gplus   - cell array, upper SRG boundary (Gauss plane)
%           .gminus  - cell array, lower SRG boundary (Gauss plane)
%           .rawdata - table from srg_compute (.freq, .maxphase, etc.)
%
%   Example - overlay three 3D SRG surfaces at τ = 1, 0.6, 0.3:
%       G  = tf([1], [1 1]);
%       hd = srg_homotopy(G, [1.0, 0.6, 0.3], -2, 3, 200, 64);
%       figure;
%       for k = 1:numel(hd)
%           srg_plot_3d(hd(k).rawdata, hd(k).gplus, 0, 0, k);
%           hold on
%       end
%       legend(arrayfun(@(s) sprintf('\\tau = %.2g', s.tau), hd, ...
%                       'UniformOutput', false));
%
%   Example - 2D comparison of homotopy slices:
%       hd = srg_homotopy(G, linspace(0.2, 1.0, 5), -2, 3, 200, 64);
%       srg_plot_compare(hd.gplus, ...
%           'Names', arrayfun(@(s) sprintf("\\tau=%.2g",s.tau), hd));
%
%   See also SRG_COMPUTE, SRG_PLOT_3D, SRG_PLOT_COMPARE.

    arguments
        G
        tau_values  (1,:) double
        rangemin    (1,1) double
        rangemax    (1,1) double
        estpoints   (1,1) double {mustBePositive, mustBeInteger}
        points      (1,1) double {mustBePositive, mustBeInteger}
        options.FreqScale string {mustBeMember(options.FreqScale,["log","linear"])} = "log"
    end

    n_tau = numel(tau_values);

    %% Preallocate struct array
    homotopy_data = struct('tau',     cell(1, n_tau), ...
                           'W',       cell(1, n_tau), ...
                           'A',       cell(1, n_tau), ...
                           'gplus',   cell(1, n_tau), ...
                           'gminus',  cell(1, n_tau), ...
                           'rawdata', cell(1, n_tau));

    %% One full freq-wise SRG per τ
    for k = 1:n_tau
        tau   = tau_values(k);
        G_tau = tau * G;

        [W, A, gp, gm, rd] = srg_compute(G_tau, rangemin, rangemax, ...
                                          estpoints, points, ...
                                          'FreqScale', options.FreqScale);

        homotopy_data(k).tau     = tau;
        homotopy_data(k).W       = W;
        homotopy_data(k).A       = A;
        homotopy_data(k).gplus   = gp;
        homotopy_data(k).gminus  = gm;
        homotopy_data(k).rawdata = rd;
    end
end
