function out = srg_qsr_multiplier(gplus, rawdata, options)
%SRG_QSR_MULTIPLIER  Frequency-wise QSR multiplier from SRG bounding circles.
%
%   OUT = SRG_QSR_MULTIPLIER(GPLUS, RAWDATA) computes, at each frequency,
%   the tightest circular multiplier Pi(c,r) containing the SRG slice,
%   then fits rational transfer functions to the multiplier parameters.
%
%   The circular multiplier region for disk D^int(c,r) is (Proposition 1,
%   arXiv:2604.05567):
%
%       Pi_int(c,r) = [-1,      c     ]    Q = -1, S = c, R = r^2-c^2
%                     [ c,  r^2 - c^2 ]
%
%   Pi_int is positive-negative (soft-hard equivalent) iff r > |c|.
%
%   Rational fitting uses Sanathanan-Koerner (SK) weighted least-squares
%   iterations.  Three weight layers are combined:
%
%     (1) SK weight        1/|A(jw)|   -- removes linearisation bias
%     (2) Frequency weight 1/w         -- equalises log-spaced data
%     (3) Relative weight  1/|y_i|     -- minimises relative error
%
%   The TF order is selected automatically by minimising AIC over 1/MaxOrder.
%
%   Inputs:
%       gplus   - Cell array of SRG upper boundaries from SRG_COMPUTE
%       rawdata - Table with .freq column (Hz) from SRG_COMPUTE
%
%   Name-Value Arguments:
%       MaxOrder - Maximum TF order to search (default: 8)
%       SKIter   - Number of SK iterations per order (default: 20)
%       Plot     - Plot per-frequency data and TF fits (default: false)
%
%   Outputs (struct):
%       out.freq     - Frequency vector (Hz)
%       out.c        - Per-frequency bounding-circle center (real)
%       out.r        - Per-frequency bounding-circle radius
%       out.R        - Per-frequency R parameter: r^2 - c^2
%       out.pos_neg  - Logical: true where Pi_int is positive-negative (r > |c|)
%       out.S_tf     - Fitted TF for c(omega)  [S parameter]
%       out.R_tf     - Fitted TF for R(omega)  [R parameter]
%       out.S_ord    - Selected order for S_tf
%       out.R_ord    - Selected order for R_tf
%       out.Pi       - Struct with fields Q=-1, S=S_tf, R=R_tf
%
%   Example:
%       G = tf([1],[1 2 1]);
%       [~,~,gp,~,rd] = srg_compute(G,-2,3,200,32);
%       out = srg_qsr_multiplier(gp, rd, 'Plot', true);
%       fprintf('Selected orders: S=%d  R=%d\n', out.S_ord, out.R_ord);
%
%   See also SRG_BOUNDING_CIRCLE, SRG_COMPUTE

    arguments
        gplus    (1,:) cell
        rawdata  table
        options.MaxOrder (1,1) double {mustBePositive, mustBeInteger} = 8
        options.SKIter   (1,1) double {mustBePositive, mustBeInteger} = 20
        options.Plot     (1,1) logical = false
    end

    % -----------------------------------------------------------------------
    %  Per-frequency bounding circles
    % -----------------------------------------------------------------------
    nf    = numel(gplus);
    c_vec = zeros(nf, 1);
    r_vec = zeros(nf, 1);

    % SISO detection: content-based (same curve repeated at every frequency)
    is_siso = (nf > 1) && ...
        (max(abs(gplus{1}(:) - gplus{2}(:))) < 1e-10 * (max(abs(gplus{1}(:))) + eps));

    if is_siso
        [c0, r0] = srg_bounding_circle(gplus{1}(:));
        c_vec(:) = c0;
        r_vec(:) = r0;
    else
        for ii = 1:nf
            pts = gplus{ii}(:);
            if numel(pts) < 2 || max(abs(pts)) < eps
                c_vec(ii) = 0;
                r_vec(ii) = 0;
            else
                [c_vec(ii), r_vec(ii)] = srg_bounding_circle(pts);
            end
        end
    end

    R_vec   = r_vec.^2 - c_vec.^2;
    pos_neg = r_vec > abs(c_vec);

    % -----------------------------------------------------------------------
    %  Fit rational TFs via SK iteration with AIC order selection
    % -----------------------------------------------------------------------
    w = 2 * pi * rawdata.freq(:);      % rad/s column vector

    [S_tf, s_ord] = sk_fit_aic(c_vec, w, options.MaxOrder, options.SKIter);
    [R_tf, r_ord] = sk_fit_aic(R_vec, w, options.MaxOrder, options.SKIter);

    % -----------------------------------------------------------------------
    %  Signed-absolute overbound: guarantee |R_tf(jw)| >= |R_i| with the
    %  same sign at every frequency.
    %
    %  Required condition:
    %      R_i > 0  =>  R_tf(jw_i) >= R_i   (fit above positive data)
    %      R_i < 0  =>  R_tf(jw_i) <= R_i   (fit below negative data)
    %
    %  A constant additive shift cannot satisfy both simultaneously:
    %  adding a positive delta corrects the positive region but pushes
    %  the negative region toward zero (wrong direction).
    %
    %  Fix: multiplicative scaling of the numerator by k >= 1.
    %      R_tf_corrected(s) = k * R_tf(s)   [multiply numerator by k]
    %  This stretches both the positive and negative parts away from zero:
    %      R_i > 0, k*R_fit > R_fit > 0  =>  larger positive  (correct)
    %      R_i < 0, k*R_fit < R_fit < 0  =>  more negative    (correct)
    %
    %  k is chosen as max_i |R_i| / |R_fit_i|, computed only over
    %  frequencies where (a) the fit has the correct sign and (b) |R_i|
    %  is above a noise floor (1% of max|R|) to avoid division near zero.
    %
    %  If the fit has the wrong sign at any frequency, MaxOrder should be
    %  increased; a warning is issued and k is computed from valid points.
    % -----------------------------------------------------------------------
    R_fit = real(squeeze(freqresp(R_tf, w)));

    tol          = 1e-2 * max(abs(R_vec));
    signs_match  = sign(R_fit) == sign(R_vec);
    above_tol    = abs(R_vec) > tol;
    valid        = signs_match & above_tol;

    n_mismatch = sum(~signs_match & above_tol);
    if n_mismatch > 0
        warning('srg_qsr_multiplier:signMismatch', ...
            ['R_tf has wrong sign at %d/%d frequencies. ' ...
             'Consider increasing MaxOrder.'], n_mismatch, nf);
    end

    if any(valid)
        ratios = abs(R_vec(valid)) ./ max(abs(R_fit(valid)), eps);
        k = max(1, max(ratios));
    else
        k = 1;
    end

    if k > 1
        [b_R, a_R] = tfdata(R_tf, 'v');
        R_tf = tf(k * b_R, a_R);
    end

    % -----------------------------------------------------------------------
    %  Pack output
    % -----------------------------------------------------------------------
    out.R_scale = k;        % multiplicative overbound factor applied
    out.freq    = rawdata.freq(:);
    out.c       = c_vec;
    out.r       = r_vec;
    out.R       = R_vec;
    out.pos_neg = pos_neg;
    out.S_tf    = S_tf;
    out.R_tf    = R_tf;
    out.S_ord   = s_ord;
    out.R_ord   = r_ord;
    out.Pi.Q    = -1;
    out.Pi.S    = S_tf;
    out.Pi.R    = R_tf;

    % -----------------------------------------------------------------------
    %  Optional diagnostic plot
    % -----------------------------------------------------------------------
    if options.Plot
        S_fit = real(squeeze(freqresp(S_tf, w)));
        R_fit = real(squeeze(freqresp(R_tf, w)));
        freq  = out.freq;

        figure;

        subplot(3,1,1);
        semilogx(freq, c_vec, 'b.', 'DisplayName', 'c_i  (data)');
        hold on;
        semilogx(freq, S_fit, 'r-', 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('S(s) fit  (order %d)', s_ord));
        grid minor;
        legend('Location', 'best');
        ylabel('Center  c(\omega)');
        title('QSR Multiplier Parameters');

        subplot(3,1,2);
        R_fit_plot = real(squeeze(freqresp(R_tf, w)));
        semilogx(freq, R_vec,      'b.', 'DisplayName', 'R_i  (data)');
        hold on;
        semilogx(freq, R_fit_plot, 'r-', 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('R(s) fit + overbound  (order %d, \\delta=%.3g)', ...
                                        r_ord, out.R_scale));
        grid minor;
        legend('Location', 'best');
        ylabel('R(\omega) = r^2 - c^2');

        subplot(3,1,3);
        semilogx(freq, r_vec,      'b-',  'LineWidth', 1.5, 'DisplayName', 'radius r');
        hold on;
        semilogx(freq, abs(c_vec), 'r--', 'LineWidth', 1.5, 'DisplayName', '|c|');
        semilogx(freq(pos_neg), r_vec(pos_neg), 'g.', 'MarkerSize', 6, ...
                 'DisplayName', 'pos-neg satisfied (r > |c|)');
        grid minor;
        legend('Location', 'best');
        ylabel('Magnitude');
        xlabel('Frequency (Hz)');
        title('Positive-negative condition  r > |c|');
    end

end

% ==========================================================================
%  sk_fit_aic  --  SK iteration with AIC-based order selection
%
%  Fits  H(jw) ~ B(jw)/A(jw)  to real-valued data y at frequencies w
%  (w in rad/s, log-spaced).
%
%  Polynomial coefficients are real; enforced by splitting the complex
%  regression into real and imaginary parts:
%
%      Re{Phi} * x  =  y       (nf equations)
%      Im{Phi} * x  =  0       (nf equations)
%
%  where  Phi = [pow_b , -diag(y)*pow_a]  and  x = [b_0...b_n, a_1...a_n].
%
%  Three weight layers are multiplied together at each SK step:
%
%    (1) SK weight        1/|A(jw_i)|
%            Removes the bias introduced by linearising B=yA.
%            This is the classical Sanathanan-Koerner correction.
%
%    (2) Frequency weight 1/w_i
%            For log-spaced data, linear LS over-represents high frequencies
%            (many points per decade at high w, few at low w).  Dividing by
%            w_i restores equal weight per log-decade, matching the visual
%            quality metric on a Bode/semilogx plot.
%
%    (3) Global amplitude scale  1/max(|y|)
%            Normalises the LS problem so the order penalty in AIC is
%            commensurable across S and R fits.  Per-point 1/|y_i| was
%            originally used for relative-error minimisation, but it blows
%            up wherever R(omega) = r^2 - c^2 crosses zero, causing the
%            SK iteration to collapse to a near-constant fit at zero and
%            AIC to reward order 1 spuriously.
%
%  AIC with relative MSE:
%      AIC = log(rel_MSE) + 2*(2n+1)/nf
%  where rel_MSE = mean( ((y-y_hat)/max(|y|,y_scale))^2 ).
%  Using relative MSE makes the order penalty commensurable across the
%  S and R fits regardless of their absolute magnitudes.
% ==========================================================================
function [sys, best_ord] = sk_fit_aic(y, w, max_ord, sk_iter)

    nf      = numel(w);
    z       = 1j * w(:);           % imaginary-axis evaluation points
    y       = y(:);                 % column, real-valued
    w       = w(:);

    % Frequency weight: 1/omega equalises log-spaced data (one weight per
    % log-decade regardless of how many points fall in that decade).
    w_freq  = 1 ./ w;

    % Global amplitude scale: normalises the LS problem without per-point
    % relative weighting.  Per-point 1/|y_i| blows up at zero crossings
    % (R = r^2 - c^2 crosses zero) and destabilises the SK iteration.
    y_global = max(abs(y)) + eps;
    w_rel    = ones(nf,1) / y_global;             % same for all points

    w_base  = w_freq .* w_rel;
    w_base  = w_base / mean(w_base);              % normalise to unit mean

    best_aic = Inf;
    best_ord = 1;
    sys      = tf(mean(y), 1);      % fallback: constant

    for n = 1:max_ord

        % Basis: B(z) = b_0 + b_1*z + ... + b_n*z^n  (monic A: a_0 = 1)
        pow_b = bsxfun(@power, z, 0:n);    % nf x (n+1)
        pow_a = bsxfun(@power, z, 1:n);    % nf x  n

        % Linearised regression matrix
        Phi    = [pow_b, -bsxfun(@times, y, pow_a)];   % nf x (2n+1)

        % Real/imaginary split (enforces real coefficients)
        Phi_ri = [real(Phi); imag(Phi)];   % 2*nf x (2n+1)
        rhs    = [y; zeros(nf,1)];         % 2*nf x 1

        % Initialise SK weights with base weights only
        wts = repmat(w_base, 2, 1);
        x   = zeros(2*n+1, 1);

        for iter = 1:sk_iter
            Phi_w = bsxfun(@times, wts, Phi_ri);
            rhs_w = wts .* rhs;

            x = Phi_w \ rhs_w;

            % Update: multiply base weight by SK correction 1/|A(jw)|
            a_tail  = x(n+2:end);
            A_resp  = 1 + pow_a * a_tail;          % nf x 1, complex
            w_sk    = 1 ./ (abs(A_resp) + eps);    % SK layer
            w_sk    = w_sk / mean(w_sk);            % keep scale stable
            w_total = w_base .* w_sk;
            wts     = repmat(w_total, 2, 1);
        end

        b_coef = real(x(1:n+1)).';           % ascending: b_0 b_1 ... b_n
        a_coef = real([1; x(n+2:end)]).';    % ascending: 1  a_1 ... a_n

        % Evaluate nonlinear residual for AIC
        B_resp = pow_b * b_coef(:);
        A_resp = 1 + pow_a * a_coef(2:end).';
        y_hat  = real(B_resp ./ A_resp);

        rel_err  = (y - y_hat) / y_global;
        rel_mse  = mean(rel_err.^2);
        aic      = log(rel_mse + eps) + 2*(2*n + 1) / nf;

        if aic < best_aic
            best_aic = aic;
            best_ord = n;
            % tf() expects descending powers: [b_n ... b_1 b_0]
            sys = tf(fliplr(b_coef), fliplr(a_coef));
        end

    end

end