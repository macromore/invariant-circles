function weightValue = weightFunction(k,numPoints)
% Compute the exponential weight value for the weighted average
%
% Define variables:
% index         - index of summand (input)
% numPoints     - number of iterates to use (input)
% unitIndex     - unitized position
% weightValue   - numerator for weight function
% 
% Dependencies: 
% N/A
    % unitize the position
    unitIndex = k/numPoints;
    % assign appropriate weight to t
    if unitIndex > 0 && unitIndex < 1
        weightValue = exp(1/(unitIndex*(unitIndex-1)));
    else % out of bounds
        weightValue = 0;
    end % end if
end % end weightFunction