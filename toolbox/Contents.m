% SRG Toolbox - Scaled Relative Graph analysis for LTI systems
% Version 1.0  2026
%
% Core computation
%   srg_compute          - Frequency-wise SRG of an LTI system or constant matrix
%   srg_hard             - Hard SRG boundary (accounts for unstable poles and NMP zeros)
%   srg_scaled_graph     - Soft SRG (Scaled Graph) via hyperbolic convex hull in BK disk
%   srg_homotopy         - Frequency-wise SRG of tau*G for multiple tau values
%
% Visualization -- Gauss plane (2-D)
%   srg_plot_gauss       - Plot the SRG in the Gauss (complex) plane
%   srg_plot_beltrami    - Plot the SRG in the Beltrami-Klein disk
%   srg_plot_compare     - Overlay multiple SRGs with automatic color cycling
%
% Visualization -- 3-D (frequency axis)
%   srg_plot_3d                  - 3D wire/filled plot of SRG boundaries vs. frequency
%   srg_plot_3d_surface          - 3D smooth surface plot of SRG boundaries
%   srg_plot_3d_compare          - Overlay two SRG surfaces with intersection highlight
%   srg_plot_3d_beltrami         - 3D wire plot in the Beltrami-Klein disk
%   srg_plot_3d_beltrami_compare - Overlay two SRG surfaces in the Beltrami-Klein disk
%   srg_plot_witness_3d          - Overlay min-distance witness segments on a 3D plot
%
% Multiplier design
%   srg_bounding_circle  - Minimum enclosing circle for an SRG boundary set
%   srg_qsr_multiplier   - Frequency-wise (Q,S,R) multiplier from SRG bounding circles
%
% Stability and margins
%   srg_stability_margin - Frequency-wise stability margin (min distance between two SRGs)
%
% Style and export
%   srg_style            - Unified style configuration struct (themes: default/dark/publication)
%   srg_apply_style      - Apply SRG style configuration to the current axes
%   srg_export           - Publication-ready export: PNG (always), TikZ (vector), PDF/EPS (best-effort)
%
% Internal utilities (bg/)
%   beltrami_map_matrix  - Beltrami-Klein mapping for a matrix operator
%   beltrami_map_scalar  - Beltrami-Klein mapping for scalar values (eigenvalues)
%   beltrami_inv         - Inverse Beltrami-Klein mapping to the Gauss plane
%   srg_to_polyshape     - Convert complex SRG boundary vector to a polyshape
%   field_of_values      - Field of values (numerical range) of a matrix
%
% Examples
%   example_01_freqwise_srg - SISO and MIMO frequency-wise SRG, stability margin, witness
%   example_02_soft_srg     - Soft SRG (hyperbolic convex hull) for SISO and MIMO
%   example_03_hard_srg     - Hard SRG for systems with unstable poles or NMP zeros
%   example_04_qsr_export   - QSR multiplier extraction and publication figure export
%
% See also SETUP
