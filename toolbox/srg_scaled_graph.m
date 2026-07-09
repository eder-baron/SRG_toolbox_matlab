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
%       sg_plus  - Upper branch of the Scaled Graph boundary (Gauss plane)
%       sg_minus - Lower branch of the Scaled Graph boundary (Gauss plane)
%       hull_bk  - Complex vector of the convex hull boundary in the
%                  Beltrami-Klein disk
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

    hull_idx = convhull(x, y);   % returns closed index (first == last)
    hull_x   = x(hull_idx);
    hull_y   = y(hull_idx);

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

    % Close the boundary
    hull_refined_x(end+1) = hull_refined_x(1);
    hull_refined_y(end+1) = hull_refined_y(1);

    hull_bk = hull_refined_x + 1i * hull_refined_y;

    %% ------------------------------------------------------------------ %
    %  Step 3: Map hull boundary to the Gauss plane via g^{-1}
    %% ------------------------------------------------------------------ %
    [sg_plus, sg_minus] = beltrami_inv(hull_bk);
    sg_p{1}=sg_plus;
    sg_m{1}=sg_minus;
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