function [ps, x_clean, y_clean] = srg_to_polyshape(srg)
%SRG_TO_POLYSHAPE  Convert complex SRG boundary to a clean polyshape.
%
%   [PS, X_CLEAN, Y_CLEAN] = SRG_TO_POLYSHAPE(SRG) takes a complex
%   vector of SRG boundary points and produces a valid MATLAB polyshape,
%   preserving the original boundary traversal order.
%
%   Outputs:
%       ps      - polyshape object
%       x_clean - cleaned real parts (original order)
%       y_clean - cleaned imaginary parts (original order)

    x = real(srg);
    y = imag(srg);

    % Remove NaN and Inf
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    % Remove consecutive near-duplicates (preserve order)
    if length(x) > 1
        tol = max(1e-15, max(sqrt(x.^2 + y.^2)) * 1e-10);
        keep = [true; sqrt(diff(x).^2 + diff(y).^2) > tol];
        x = x(keep);
        y = y(keep);
    end

    x_clean = x(:);
    y_clean = y(:);

    if length(x_clean) < 3
        ps = polyshape();
        return;
    end

    % Close if not already closed
    if sqrt((x_clean(end)-x_clean(1))^2 + (y_clean(end)-y_clean(1))^2) > tol
        x_clean(end+1) = x_clean(1);
        y_clean(end+1) = y_clean(1);
    end

    % Create polyshape preserving original order
    ps = polyshape(x_clean, y_clean, 'Simplify', false);

end
