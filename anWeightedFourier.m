function [an , bn] = anWeightedFourier(orbit,rho,n,totalWeight)
% Compute the n-th complex Fourier coefficient of a trajectory using
% weighted Birkhoff ergodic theorem. 
% 
% Define variables:
% an          - partial sum for the n-th coefficient in the x-series
% bn          - partial sum for the n-th coefficient in the y-series
% n           - index of the desired coeffieients (input)
% orbit       - initial data, the trajectory (input)
% rho         - computed rotation number (input)
% totalWeight - denominator for weight function (input)
% weightValue - numerator for weight function
% 
% Dependencies: 
% weightFunction.m
% Initialize partial sums
    an = 0; 
    bn = 0;
    % Compute the sum portion of the weighted Birkhoff averages
    for j = 1:length(orbit)-1
        weightValue = weightFunction(j,length(orbit)-1);
        an = an + orbit(1,j+1)*exp(1i*2*pi*j*n*rho)*weightValue;
        bn = bn + orbit(2,j+1)*exp(1i*2*pi*j*n*rho)*weightValue;
    end % end for loop
    % Divide by the total weight
    an = an/totalWeight;
    bn = bn/totalWeight;
end % end function