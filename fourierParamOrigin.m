function [a, b] = fourierParamOrigin(orbit, N,K,rho)
% A possibly defunct script to compute parameterizations of a trajectory,
% period one only when origin is inside the orbit
% 
% Define variables:
% orbit         - initial data, the trajectory (input)
% N             - number of Fourier modes (input)
% K             - number of closed orbits (input)
% rho           - computed rotation number (input)
% fouriercoeffs - array to store coefficients
% a             - x coordinate fourier series
% b             - y coordinate fourier series
% 
% Dependencies: 
% Fourier.m
% anWeightedFourier.m
% 
% Notes: Doesn't need to involve K
    % set default value of N and K
    if nargin == 1
        N = 20;
        K = 1;
    elseif nargin == 2
        K = 1;
    end % end if
    % Correct rho for higher period
    rho = rho*K;
    % iniitialize coefficients array
    fouriercoeffs = zeros(2,2*N+1);
    % Compute Total Weight
    totalWeight = 0;
    for k = 1:length(orbit)
        totalWeight = totalWeight + weightFunction(k,length(orbit)-1);
    end % end for loop
    % Compute coefficients
    for n = -N:N
        %[an, bn] = anFourier(orbit,rho,n);
        [an, bn] = anWeightedFourier(orbit,rho,n,totalWeight);
        fouriercoeffs(:,n+N+1) = [an ; bn];
    end % end for loop
    % Convert to Fourier class
    a = Fourier(fouriercoeffs(1,:));
    b = Fourier(fouriercoeffs(2,:));
end % end fourierParamOrigin