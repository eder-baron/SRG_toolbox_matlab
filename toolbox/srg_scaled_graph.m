function [sg_p, sg_m, hull_bk] = srg_scaled_graph(input, options)
%SRG_SCALED_GRAPH  Scaled Graph via hyperbolic convex hull.
%
%   [SG_PLUS, SG_MINUS, HULL_BK] = SRG_SCALED_GRAPH(W) takes the
%   frequency-wise Beltrami-Klein numerical range boundaries from
%   srg_compute and returns the Scaled Graph: the hyperbolic convex hull
%   of all frequency-wise Scaled Relative Graphs.
%
%   [SG_PLUS, SG_MINUS, HULL_BK] = SRG_SCALED_GRAPH(GPLUS, GMINUS=GMINUS)
%   accepts the SRG directly in the Gauss plane (gplus/gminus cell arrays
%   from srg_compute), maps them forward into the BK disk, and then
%   computes the Scaled Graph. Use this form when you already have
%   gplus/gminus and do not need to recompute W.
%
%   The algorithm:
%       BK-disk input path:
%           1. Pools all BK-disk boundary points across frequencies
%           2. Computes their Euclidean convex hull in the BK disk
%           3. Refines the hull boundary by edge interpolation
%           4. Maps the hull boundary to the Gauss plane via g^{-1}
%
%       Gauss-plane input path:
%           1. Maps gplus and gminus forward into the BK disk via the
%              scalar Beltrami-Klein map  f(lambda) = (conj(l)-i)(l-i)
%                                                      ----------------
%                                                        1 + conj(l)*l
%           2-4. Same as above.
%
%   Since the Beltrami-Klein model represents hyperbolic geometry with
%   Euclidean straight lines as geodesics, the Euclidean convex hull in
%   the BK disk IS the hyperbolic convex hull.
%
%   Degenerate (collinear) point sets: for some systems the pooled
%   BK-disk points are exactly collinear -- e.g. a first-order SISO
%   system, whose Nyquist curve's BK image is a straight chord. A
%   proper 2D convex hull is undefined for collinear input (CONVHULL
%   errors), so this case is detected up front via the rank of the
%   centered point set, and the two extreme points along the principal
%   direction are used directly as a degenerate hull: the single open
%   chord between them (not closed back to the start, since a line
%   segment has no interior to enclose).
%
%   Inputs (positional):
%       W      - Cell array of complex vectors, each containing BK-disk
%                boundary points at one frequency (output of srg_compute).
%                Used when 'InputDomain' is 'BK' (default).
%
%   Name-Value Arguments:
%       InputDomain      - 'BK' (default) or 'Gauss'.
%                          When 'Gauss', supply GPlus and GMinus as well.
%       GPlus            - Cell array of gplus curves (Gauss plane).
%                          Required when InputDomain = 'Gauss'.
%       GMinus           - Cell array of gminus curves (Gauss plane).
%                          Required when InputDomain = 'Gauss'.
%       RefinementPoints - Number of interpolation points per hull edge
%                          (default: 100).
%
%   Outputs:
%       sg_plus  - Upper branch of the Scaled Graph boundary (Gauss plane),
%                  wrapped in a 1x1 cell for direct compatibility with the
%                  Gauss-plane plot functions (srg_plot_gauss etc.)
%       sg_minus - Lower branch of the Scaled Graph boundary (Gauss plane),
%                  wrapped in a 1x1 cell (see sg_plus)
%       hull_bk  - 1x1 cell containing the complex vector of the convex
%                  hull boundary in the Beltrami-Klein disk. Wrapped in a
%                  cell (rather than returned as a bare vector) so it is
%                  directly compatible with srg_plot_beltrami and
%                  srg_plot_3d_beltrami, both of which expect a cell
%                  array of BK-plane curves and already have a
%                  dedicated code path for a single static curve
%                  (numel(W)==1): srg_plot_beltrami dispatches on
%                  isa(TransferFcn,'double') for that path, so pass any
%                  double as TransferFcn (its value is unused there);
%                  srg_plot_3d_beltrami dispatches on numel(gplus1)==1
%                  directly, so no such placeholder is needed. See the
%                  examples below.
%
%   Examples:
%       % --- From BK-disk output of srg_compute ---
%       G = tf([1 2], [1 3 2]);
%       [W, ~, gplus, gminus, ~] = srg_compute(G, -2, 3, 200, 32);
%       [sg_plus, sg_minus, hull_bk] = srg_scaled_graph(W);
%
%       % --- From Gauss-plane SRG directly ---
%       [sg_plus, sg_minus, hull_bk] = srg_scaled_graph([], ...
%           'InputDomain', 'Gauss', 'GPlus', gplus, 'GMinus', gminus);
%
%       % --- Degenerate/collinear case (e.g. G = 1/(s+1)) ---
%       G = tf(1, [1 1]);
%       [W, ~, ~, ~, ~] = srg_compute(G, -2, 3, 200, 1);
%       [sg_plus, sg_minus, hull_bk] = srg_scaled_graph(W);  % no longer errors
%
%       % --- Plotting the BK-disk hull directly ---
%       % hull_bk is already a 1x1 cell, so it drops straight into the
%       % Beltrami plotting functions with no reshaping needed.
%       figure;
%       srg_plot_beltrami(1, hull_bk, {[]}, 1, 1);   % TransferFcn=1 is a
%                                                     % placeholder double,
%                                                     % only its class matters
%       figure;
%       srg_plot_3d_beltrami(rawdata_single_freq, hull_bk, 0, 0, 1);
%
%   See also SRG_COMPUTE, SRG_PLOT_GAUSS, SRG_PLOT_BELTRAMI

    arguments
        input                % W cell array (BK) or [] when InputDomain='Gauss'
        options.InputDomain  (1,:) char  {mustBeMember(options.InputDomain, {'BK','Gauss'})} = 'BK'
        options.GPlus        (1,:) cell  = {}
        options.GMinus       (1,:) cell  = {}
        options.RefinementPoints (1,1) double {mustBePositive, mustBeInteger} = 100
    end

    %% ------------------------------------------------------------------ %
    %  Step 0: Resolve input domain and collect BK-disk points
    %% ------------------------------------------------------------------ %
    switch options.InputDomain

        case 'BK'
            % Input is already the W cell array from srg_compute
            if isempty(input) || ~iscell(input)
                error('srg_scaled_graph:badInput', ...
                    'When InputDomain=''BK'', the first argument must be a non-empty cell array W.');
            end
            all_bk = vertcat(input{:});

        case 'Gauss'
            % Map gplus / gminus forward into the BK disk
            if isempty(options.GPlus) || isempty(options.GMinus)
                error('srg_scaled_graph:missingGauss', ...
                    'When InputDomain=''Gauss'', GPlus and GMinus must be provided.');
            end
            gp_bk = cellfun(@gauss_to_bk, options.GPlus,  'UniformOutput', false);
            gm_bk = cellfun(@gauss_to_bk, options.GMinus, 'UniformOutput', false);
            all_bk = vertcat(gp_bk{:}, gm_bk{:});
    end

    % Keep only finite points strictly inside (or on) the unit disk
    all_bk = all_bk(isfinite(all_bk) & abs(all_bk) <= 1 + 1e-9);
    if numel(all_bk) < 3
        error('srg_scaled_graph:tooFewPoints', ...
            'Fewer than 3 valid BK-disk points; cannot compute convex hull.');
    end

    %% ------------------------------------------------------------------ %
    %  Step 1: Euclidean convex hull in the BK disk
    %% ------------------------------------------------------------------ %
    x = real(all_bk);
    y = imag(all_bk);

    % Detect a degenerate (collinear / rank-1) point set via the SVD of
    % the centered coordinates: a proper 2D hull needs two non-negligible
    % singular values. This is exactly the situation for e.g. a
    % first-order SISO system, whose BK-disk image is a straight chord;
    % CONVHULL cannot handle purely collinear input.
    centroid = [mean(x), mean(y)];
    centered = [x - centroid(1), y - centroid(2)];
    [~, S, V] = svd(centered, 'econ');
    sv = diag(S);
    is_degenerate = (numel(sv) < 2) || (sv(2) <= 1e-9 * max(sv(1), eps));

    if is_degenerate
        if sv(1) <= 1e-9 * max(abs([x; y]), [], 'all')
            % All pooled points coincide (zero-length hull); nothing
            % meaningful to build a Scaled Graph from.
            error('srg_scaled_graph:degenerateHull', ...
                'All BK-disk points coincide; cannot form a Scaled Graph.');
        end
        % Collinear points: take the two extreme points along the
        % principal direction and treat the segment between them as a
        % single open chord. This is a genuine 1-D object (a line has no
        % interior), so unlike CONVHULL's closed polygon output it is
        % NOT traversed back to the start -- doing so would retrace the
        % same straight line on top of itself, which renders as an
        % overlapping/doubled stroke rather than one clean line.
        proj = centered * V(:,1);
        [~, i_min] = min(proj);
        [~, i_max] = max(proj);
        hull_idx = [i_min; i_max];
        hull_x   = x(hull_idx);
        hull_y   = y(hull_idx);
    else
        hull_idx = convhull(x, y);   % returns closed index (first == last)
        hull_x   = x(hull_idx);
        hull_y   = y(hull_idx);
    end

    %% ------------------------------------------------------------------ %
    %  Step 2: Refine hull edges
    %  The inverse BK mapping g^{-1} is nonlinear, so straight edges in
    %  the BK disk become curves in the Gauss plane. Interpolating each
    %  edge densely before mapping ensures a smooth output boundary.
    %% ------------------------------------------------------------------ %
    n_refine = options.RefinementPoints;
    n_edges  = length(hull_idx) - 1;

    hull_refined_x = zeros(n_edges * n_refine, 1);
    hull_refined_y = zeros(n_edges * n_refine, 1);

    for ii = 1:n_edges
        idx_start = (ii-1) * n_refine + 1;
        idx_end   = ii     * n_refine;
        hull_refined_x(idx_start:idx_end) = linspace(hull_x(ii), hull_x(ii+1), n_refine);
        hull_refined_y(idx_start:idx_end) = linspace(hull_y(ii), hull_y(ii+1), n_refine);
    end

    % Close the boundary. Only meaningful for the non-degenerate polygon
    % hull (whose CONVHULL-derived edges already trace back to the start
    % on the final edge, so this just adds an explicit duplicate final
    % point). The degenerate open chord is a line segment, not a
    % polygon, and must NOT be closed -- see the comment above.
    if ~is_degenerate
        hull_refined_x(end+1) = hull_refined_x(1);
        hull_refined_y(end+1) = hull_refined_y(1);
    end

    hull_bk = hull_refined_x + 1i * hull_refined_y;

    %% ------------------------------------------------------------------ %
    %  Step 3: Map hull boundary to the Gauss plane via g^{-1}
    %% ------------------------------------------------------------------ %
    [sg_plus, sg_minus] = beltrami_inv(hull_bk);
    sg_p{1}=sg_plus;
    sg_m{1}=sg_minus;

    % Wrap in a 1x1 cell to match sg_p/sg_m's convention and to plug
    % directly into srg_plot_beltrami / srg_plot_3d_beltrami, which both
    % expect a cell array of BK-plane curves (see the Outputs doc above).
    hull_bk = {hull_bk};
end

%% ========================================================================
%  Local helper: forward Beltrami-Klein map  f: Gauss plane -> BK disk
%  For a scalar complex value lambda (a point in the Gauss plane / C),
%  the map is the scalar version of ftop / f used in srg_compute:
%
%      f(lambda) = (conj(lambda) - i)(lambda - i)
%                  --------------------------------
%                        1 + conj(lambda)*lambda
%
%  This is exactly the scalar case of ftop.m and f.m in the original code.
%  Points at lambda = i (the "north pole") map to the boundary |w|=1 and
%  are clipped.
%% ========================================================================
function w = gauss_to_bk(lambda)
%GAUSS_TO_BK  Map Gauss-plane SRG points into the Beltrami-Klein disk.
%   Lambda may be a vector of complex numbers.
    lambda  = lambda(:);
    cl      = conj(lambda);
    denom   = 1 + cl .* lambda;        % = 1 + |lambda|^2  (real, >= 1)
    w       = (cl - 1i) .* (lambda - 1i) ./ denom;
    % Clip any numerically out-of-range values (should not occur for
    % well-formed SRG points, but guard against floating-point overshoot)
    w(abs(w) > 1) = w(abs(w) > 1) ./ abs(w(abs(w) > 1));
end
