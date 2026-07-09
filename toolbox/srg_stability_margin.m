function [stab_margin, freq, x1, x2] = srg_stability_margin(gplus1, gplus2, rawdata1)
%SRG_STABILITY_MARGIN  Frequency-wise stability margin between two SRGs.
%
%   [STAB_MARGIN, FREQ, X1, X2] = SRG_STABILITY_MARGIN(gplus1, gplus2, rawdata1)
%   computes the minimum distance between two SRG boundaries at each
%   frequency point, providing a frequency-resolved stability margin.
%
%   Inputs:
%       gplus1   - Cell array of SRG boundary for system 1
%       gplus2   - Cell array of SRG boundary for system 2
%       rawdata1 - Frequency vector (or table .freq column)
%
%   Outputs:
%       stab_margin - Minimum distance at each frequency
%       freq        - Frequency vector
%       x1          - Closest point on SRG 1 at each frequency
%       x2          - Closest point on SRG 2 at each frequency
%
%   See also SRG_COMPUTE,

    freq = rawdata1;
    a = length(freq);
    stab_margin = zeros(a, 1);
    x1 = NaN(a, 1);
    x2 = NaN(a, 1);

    is_scalar1 = isscalar(gplus1);
    is_scalar2 = isscalar(gplus2);

    for ii = 1:a
        % Select appropriate cell index
        idx1 = min(ii, length(gplus1));
        idx2 = min(ii, length(gplus2));
        if is_scalar1, idx1 = 1; end
        if is_scalar2, idx2 = 1; end

        re1 = real(gplus1{idx1}(:));  im1 = imag(gplus1{idx1}(:));
        re2 = real(gplus2{idx2}(:));  im2 = imag(gplus2{idx2}(:));

        % Pairwise Euclidean distances 
        dre = bsxfun(@minus, re1, re2.');
        dim = bsxfun(@minus, im1, im2.');
        allDistances = sqrt(dre.^2 + dim.^2);

        [stab_margin(ii), index] = min(allDistances(:));

        if isnan(stab_margin(ii))
            stab_margin(ii) = 100;
        end

        [x1_index, x2_index] = ind2sub(size(allDistances), index);
        x1(ii) = gplus1{idx1}(x1_index);
        x2(ii) = gplus2{idx2}(x2_index);
    end

end