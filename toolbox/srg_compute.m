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
%   INVERSE operator, SRG(T^-1), computed by taking the pointwise
%   inverse 1/z of every point z on the Gauss-plane SRG boundary
%   (gplus, gminus). This uses the identity SRG(T^-1) = {1/z : z in
%   SRG(T)}, which holds whenever 0 is not in the interior of SRG(T).
%   The inverted Gauss-plane points are then mapped back into the
%   Beltrami-Klein disk via the forward BK map, and returned as W, so
%   that W also reflects the inverse operator rather than the
%   original one. A (the eigenvalues in the BK disk) is left
%   unchanged, since it is not derived from gplus/gminus.
%
%   At every frequency slice, srg_compute checks whether the origin
%   lies in the interior of the ORIGINAL (non-inverted) SRG. If it
%   does, the pointwise-inverted region is unbounded and contains the
%   point at infinity; a warning is issued listing the frequency
%   range(s) where this occurs, and rawdata.zero_interior flags the
%   affected frequency indices.
%
%   Inputs:
%       TransferFcn - Transfer function (tf/ss object) or constant matrix
%       rangemin    - Minimum frequency exponent (log) or value (linear)
%       rangemax    - Maximum frequency exponent (log) or value (linear)
%       estpoints   - Number of frequency evaluation points
%       points      - Angular resolution for field of values computation.
%                     Increase (e.g. 64, 128) for smoother boundaries on
%                     systems with sharply curved numerical ranges.
%
%   Name-Value Arguments:
%       FreqScale   - Frequency spacing: "log" (default) or "linear"
%       Inverse     - If true, return the pointwise-inverted (1/z) SRG,
%                     i.e. SRG(T^-1). Default false.
%
%   Outputs:
%       W       - Cell array of numerical range boundaries (Beltrami-Klein).
%                 When Inverse=true, contains the BK-disk image of the
%                 pointwise-inverted SRG instead of the original operator's
%                 numerical range.
%       A       - Cell array of eigenvalues in Beltrami-Klein plane
%       gplus   - Cell array of SRG upper boundaries (Gauss plane)
%       gminus  - Cell array of SRG lower boundaries (Gauss plane)
%       rawdata - Table with freq, maxphase, minphase, norminf, norminfm.
%                 When Inverse=true, also includes zero_interior, a
%                 logical column flagging frequencies at which the
%                 origin was interior to the ORIGINAL SRG (i.e. where
%                 the inverse SRG is unbounded and contains infinity).
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
%       [W,A,gp,gm,rd] = srg_compute(G, -2, 3, 200, 32, 'Inverse', true);
%       if any(rd.zero_interior)
%           disp('Inverse SRG is unbounded (contains infinity) at some frequencies.')
%       end
%
%   See also SRG_SCALED_GRAPH, SRG_PLOT_GAUSS, SRG_PLOT_BELTRAMI, SRG_HOMOTOPY

    arguments
        TransferFcn
        rangemin    (1,1) double = NaN
        rangemax    (1,1) double = NaN
        estpoints   (1,1) double {mustBePositive, mustBeInteger} = 200
        points      (1,1) double {mustBePositive, mustBeInteger} = 64
        options.FreqScale string {mustBeMember(options.FreqScale, ...
            ["log","linear","auto"])} = "auto"
        options.Inverse (1,1) logical = false
    end

    %% Frequency vector
    auto = options.FreqScale == "auto" || isnan(rangemin) || isnan(rangemax);

    if auto && ~isa(TransferFcn, 'double')
        [~, ~, wout] = nyquist(TransferFcn);   % rad/s, no plot
        y1 = wout(:)' / (2*pi);                % convert to Hz
        estpoints = numel(y1);
    elseif options.FreqScale == "linear"
        y1 = linspace(rangemin, rangemax, estpoints);
    else  % "log" or auto fallback for constant matrices
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

    %% Handle constant matrix: compute once, replicate across all frequencies
    if isa(TransferFcn, 'double')
        T_jw = TransferFcn;
        map_1 = beltrami_map_matrix(T_jw);
        [W1, A1] = field_of_values(map_1, 1, points, 1);
        [gp1, gm1] = beltrami_inv(W1);

        is_interior_1 = false;
        if options.Inverse
            is_interior_1 = origin_in_srg(gp1, gm1);
            gp1 = 1 ./ gp1;
            gm1 = 1 ./ gm1;

            % Map the pointwise-inverted Gauss-plane boundary back into
            % the Beltrami-Klein disk, so W reflects T^-1 too.
            W1 = [beltrami_map_scalar(gp1(:)); flipud(beltrami_map_scalar(gm1(:)))];
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
            zero_interior(ii) = is_interior_1;
        end

    else
    %% Main computation loop for transfer functions
    for ii = 1:estpoints

        T_jw = freqresp(TransferFcn, y1(ii), 'Hz');

        % Beltrami-Klein mapping of the linear operator
        map_ii = beltrami_map_matrix(T_jw);%* exp(-1j * y1(ii))

        % Numerical range (field of values) in the Beltrami-Klein plane
        [W1, A1] = field_of_values(map_ii, 1, points, 1);
        W{ii} = W1;
        A{ii} = A1;

        % Inverse Beltrami-Klein mapping to Gauss plane
        [gp_ii, gm_ii] = beltrami_inv(W{ii});

        % Optionally take the pointwise inverse: SRG(T^-1) = {1/z : z in SRG(T)}
        if options.Inverse
            zero_interior(ii) = origin_in_srg(gp_ii, gm_ii);
            gp_ii = 1 ./ gp_ii;
            gm_ii = 1 ./ gm_ii;

            % Map the pointwise-inverted Gauss-plane boundary back into
            % the Beltrami-Klein disk, so W reflects T^-1 too.
            W{ii} = [beltrami_map_scalar(gp_ii(:)); flipud(beltrami_map_scalar(gm_ii(:)))];
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
                ['The origin lies in the interior of the original SRG at ' ...
                 '%d %s. The inverse SRG therefore contains the exterior ' ...
                 'of the plotted boundary and includes the point at ' ...
                 'infinity over the frequency range(s): %s.'], ...
                 numel(idx), freqWord, rangeList);
        end
    end

end

%% Local helper: is the origin interior to the SRG boundary (gplus/gminus)?
function tf = origin_in_srg(gp, gm)
    gp = gp(:);
    gm = gm(:);
    % Close the boundary loop: upper boundary forward, lower boundary reversed
    xs = [real(gp); flipud(real(gm))];
    ys = [imag(gp); flipud(imag(gm))];

    % Degenerate (near single-point) boundaries cannot enclose the origin
    if numel(xs) < 3 || (range(xs) < eps(1) && range(ys) < eps(1))
        tf = false;
        return;
    end

    tf = inpolygon(0, 0, xs, ys);
end