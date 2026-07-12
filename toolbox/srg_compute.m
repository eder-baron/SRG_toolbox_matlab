function [W,A,gplus,gminus,rawdata] = srg_compute(TransferFcn, rangemin, rangemax, estpoints, points, options)
%SRG_COMPUTE Compute the Scaled Relative Graph of a linear operator.
%
%   [W,A,gplus,gminus,rawdata] = SRG_COMPUTE(TransferFcn, rangemin,
%   rangemax, estpoints, points) computes the SRG using logarithmically
%   spaced frequencies (default).
%
%   [...] = SRG_COMPUTE(..., 'FreqScale', 'linear') uses linearly spaced
%   frequencies instead of logarithmic spacing.
%
%   [...] = SRG_COMPUTE(..., 'Inverse', true) returns the SRG of the
%   INVERSE operator, SRG(T^-1). At each frequency, the actual matrix
%   inverse of the frequency response, inv(T_jw), is computed and run
%   through the same Beltrami-Klein mapping / field-of-values / inverse
%   BK mapping pipeline used for the non-inverse case..
%
%   At every frequency slice where Inverse=true, srg_compute checks
%   whether the resulting boundary (gplus/gminus) contains any
%   non-finite points. A warning is issued listing the
%   frequency range(s) where this occurs, and rawdata.zero_interior
%   flags the affected frequency indices.
%
%
%   Inputs:
%       TransferFcn - Transfer function (tf/ss object) or constant matrix
%       rangemin    - Minimum frequency exponent (log) or value (linear)
%       rangemax    - Maximum frequency exponent (log) or value (linear)
%       estpoints   - Number of frequency evaluation points
%       points      - Angular resolution for field of values computation.
%                     Increase (e.g. 64, 128) for smoother boundaries on
%                     systems with sharply curved numerical ranges. Also
%                     used as the chord-refinement density (see
%                     'Refine' below).
%
%   Name-Value Arguments:
%       FreqScale   - Frequency spacing: "log" (default) or "linear"
%       Inverse     - If true, return SRG(T^-1), computed directly from
%                     inv(T_jw) at each frequency. Default false.
%                     Requires TransferFcn to be square and invertible
%                     (in the freqresp sense) across the frequency
%                     sweep; behaves the same as manually computing
%                     inv(TransferFcn) and passing that in.
%       Refine      - If true (default), subdivide long BK-plane chords
%                     via SRG_REFINE_CHORDS before mapping to the Gauss
%                     plane. See above.
%       RefineTol   - Maximum real- or imaginary-part step between
%                     consecutive BK-plane boundary points before
%                     SRG_REFINE_CHORDS subdivides that chord. Default
%                     0.05. Only used when Refine=true.
%
%   Outputs:
%       W       - Cell array of numerical range boundaries (Beltrami-Klein).
%                 When Inverse=true, contains the BK-disk numerical
%                 range of inv(T_jw) instead of the original operator's.
%       A       - Cell array of eigenvalues in Beltrami-Klein plane.
%                 When Inverse=true, contains the eigenvalues of the
%                 mapped inv(T_jw) instead of the original operator's.
%       gplus   - Cell array of SRG upper boundaries (Gauss plane)
%       gminus  - Cell array of SRG lower boundaries (Gauss plane)
%       rawdata - Table with freq, maxphase, minphase, norminf, norminfm.
%                 When Inverse=true, also includes zero_interior, a
%                 logical column flagging frequencies at which the
%                 inverse SRG boundary contains non-finite (unbounded)
%                 points.
%
%   Examples:
%       % Logarithmic spacing (default)
%       G = tf([1],[1 1]);
%       [W,A,gp,gm,rd] = srg_compute(G, -2, 3, 200, 32);
%
%       % Linear spacing
%       [W,A,gp,gm,rd] = srg_compute(G, 0.01, 100, 200, 32, ...
%                                      'FreqScale', 'linear');
%
%       % Inverse SRG, SRG(G^-1)
%       [Wi,Ai,gpi,gmi,rdi] = srg_compute(G, -2, 3, 200, 32, 'Inverse', true);
%       if any(rdi.zero_interior)
%           disp('Inverse SRG is unbounded (contains infinity) at some frequencies.')
%       end
%
%       % Normal operator with real eigenvalues (needs chord refinement,
%       % which is on by default -- shown here just for illustration)
%       H = diag([3, -2/3, 1]);
%       [W,A,gp,gm,rd] = srg_compute(H, [], [], 1, 32, 'RefineTol', 0.02);
%
%   See also SRG_SCALED_GRAPH, SRG_PLOT_GAUSS, SRG_PLOT_BELTRAMI, ...
%            SRG_HOMOTOPY, SRG_REFINE_CHORDS

    arguments
        TransferFcn
        rangemin    (1,1) double = NaN
        rangemax    (1,1) double = NaN
        estpoints   (1,1) double {mustBePositive, mustBeInteger} = 200
        points      (1,1) double {mustBePositive, mustBeInteger} = 64
        options.FreqScale string {mustBeMember(options.FreqScale, ...
            ["log","linear","auto"])} = "auto"
        options.Inverse (1,1) logical = false
        options.Refine (1,1) logical = true
        options.RefineTol (1,1) double {mustBePositive} = 0.05
    end

    %% Frequency vector
    % 'auto' (nyquist-driven) selection only applies when the range was
    % genuinely left unspecified (rangemin/rangemax are NaN). If the
    % caller supplies an explicit range, it is honored even when
    % FreqScale is left at its default "auto" value -- otherwise
    % explicit rangemin/rangemax/estpoints would be silently discarded
    % whenever 'FreqScale' isn't also passed, and estpoints in the
    % caller's workspace would no longer match numel(W).
    rangeGiven = ~isnan(rangemin) && ~isnan(rangemax);
    auto = options.FreqScale == "auto" && ~rangeGiven;

    if auto && ~isa(TransferFcn, 'double')
        [~, ~, wout] = nyquist(TransferFcn);   % rad/s, no plot
        y1 = wout(:)' / (2*pi);                % convert to Hz
        estpoints = numel(y1);
    elseif options.FreqScale == "linear"
        y1 = linspace(rangemin, rangemax, estpoints);
    else  % "log" (explicit or auto fallback for constant matrices / unspecified range)
        if isnan(rangemin), rangemin = -2; end
        if isnan(rangemax), rangemax =  3; end
        y1 = logspace(rangemin, rangemax, estpoints);
    end



    %% Constant matrices use the full frequency vector (same SRG at every freq)

    %% Preallocate
    W      = cell(1, estpoints);
    A      = cell(1, estpoints);
    gplus  = cell(1, estpoints);
    gminus = cell(1, estpoints);

    maxphase = zeros(estpoints, 1);
    minphase = zeros(estpoints, 1);
    norminf  = zeros(estpoints, 1);
    norminfm = zeros(estpoints, 1);
    zero_interior = false(estpoints, 1);

    % Temporarily silence MATLAB's built-in singular-matrix warning when
    % inverting T_jw: a singular/near-singular T_jw is an expected and
    % meaningful situation here (it is exactly what produces an
    % unbounded inverse SRG), and is already reported via
    % rawdata.zero_interior and the consolidated warning below.
    if options.Inverse
        singWarnState = warning('off', 'MATLAB:singularMatrix');
        nearSingWarnState = warning('off', 'MATLAB:nearlySingularMatrix');
        cleanupWarn1 = onCleanup(@() warning(singWarnState)); %#ok<NASGU>
        cleanupWarn2 = onCleanup(@() warning(nearSingWarnState)); %#ok<NASGU>
    end

    %% Handle constant matrix: compute once, replicate across all frequencies
    if isa(TransferFcn, 'double')
        T_jw = TransferFcn;
        if options.Inverse
            T_jw = inv(T_jw);
        end
        map_1 = beltrami_map_matrix(T_jw);
        [W1, A1] = field_of_values(map_1, 1, points, 1);
        if options.Refine
            W1 = srg_refine_chords(W1, points, options.RefineTol);
        end
        [gp1, gm1] = beltrami_inv(W1);

        is_unbounded_1 = false;
        if options.Inverse
            is_unbounded_1 = any(~isfinite(gp1)) || any(~isfinite(gm1));
        end

        mp = max(abs(angle(gm1))) * 180/pi;
        mnp = min(abs(angle(gm1))) * 180/pi;
        ni = norm(gp1, Inf);
        nim = norm(gp1, -Inf);

        for ii = 1:estpoints
            W{ii}      = W1;
            A{ii}      = A1;
            gplus{ii}  = gp1;
            gminus{ii} = gm1;
            maxphase(ii) = mp;
            minphase(ii) = mnp;
            norminf(ii)  = ni;
            norminfm(ii) = nim;
            zero_interior(ii) = is_unbounded_1;
        end

    else
    %% Main computation loop for transfer functions
    for ii = 1:estpoints

        T_jw = freqresp(TransferFcn, y1(ii), 'Hz');

        % Optionally work with the actual matrix inverse of the
        % frequency response, so everything downstream (BK mapping,
        % field of values, inverse BK mapping) is computed exactly for
        % the inverse operator, at full angular resolution.
        if options.Inverse
            T_jw = inv(T_jw);
        end

        % Beltrami-Klein mapping of the linear operator
        map_ii = beltrami_map_matrix(T_jw);%* exp(-1j * y1(ii))

        % Numerical range (field of values) in the Beltrami-Klein plane
        [W1, A1] = field_of_values(map_ii, 1, points, 1);

        % Subdivide long BK-plane chords before the nonlinear mapping to
        % the Gauss plane. This is a no-op (returns W1 unchanged) unless
        % some consecutive boundary points are farther apart than
        % RefineTol, which mainly happens for normal operators (e.g.
        % real eigenvalues at this frequency slice).
        if options.Refine
            W1 = srg_refine_chords(W1, points, options.RefineTol);
        end

        W{ii} = W1;
        A{ii} = A1;

        % Inverse Beltrami-Klein mapping to Gauss plane
        [gp_ii, gm_ii] = beltrami_inv(W{ii});

        if options.Inverse
            % The inverse BK mapping has a singularity precisely where
            % SRG(T) comes close to containing the origin; flag any
            % non-finite boundary points that result.
            zero_interior(ii) = any(~isfinite(gp_ii)) || any(~isfinite(gm_ii));
        end

        gplus{ii}  = gp_ii;
        gminus{ii} = gm_ii;

        % Store scalar metrics for the generalized Bode plot
        maxphase(ii) = max(abs(angle(gminus{ii}))) * 180/pi;
        minphase(ii) = min(abs(angle(gminus{ii}))) * 180/pi;
        norminf(ii)  = norm(gplus{ii}, Inf);
        norminfm(ii) = norm(gplus{ii}, -Inf);
    end
    end

    %% Build output table
    rawdata = table;
    rawdata.freq     = y1';
    rawdata.maxphase = maxphase;
    rawdata.minphase = minphase;
    rawdata.norminf  = norminf;
    rawdata.norminfm = norminfm;

    if options.Inverse
        rawdata.zero_interior = zero_interior;

        if any(zero_interior)
            idx   = find(zero_interior);
            gaps  = find(diff(idx) > 1);
            starts = [idx(1); idx(gaps+1)];
            ends   = [idx(gaps); idx(end)];

            rangeStrs = strings(numel(starts), 1);
            for kk = 1:numel(starts)
                rangeStrs(kk) = sprintf('[%.6g, %.6g] Hz', ...
                    y1(starts(kk)), y1(ends(kk)));
            end
            rangeList = strjoin(rangeStrs, ', ');

            if numel(idx) == 1
                freqWord = 'frequency';
            else
                freqWord = 'frequencies';
            end

            warning('srg_compute:UnboundedInverseSRG', ...
                ['The inverse SRG boundary contains non-finite points at ' ...
                 '%d %s (T_jw is singular or near-singular there). The ' ...
                 'inverse SRG is unbounded and includes the point at ' ...
                 'infinity over the frequency range(s): %s.'], ...
                 numel(idx), freqWord, rangeList);
        end
    end

end
