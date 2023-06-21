function [a, b] = fourierParam(torus, modes, K, rho, perPt)
% A possibly defunct script to compute parameterizations of a trajectory,
% period one only.
% 
% Define variables:
% orbit       - initial data, the trajectory (input)
% N           - number of Fourier modes (input)
% K           - number of closed orbits (input)
% rho         - computed rotation number (input)
% avex        - average of the x trajectory
% avey        - average of the y trajectory
% oprime      - orbit re-centered at average
% a           - x coordinate fourier series
% b           - y coordinate fourier series
% 
% Dependencies: 
% Fourier.m
% fourierParamOrigin.m
    % Translate the orbit to enclose origin
    tprime = [torus(1,:)-perPt(1) ; torus(2,:) - perPt(2)];
    % Send off the translated orbit
    [a, b] = fourierParamOrigin(tprime, modes, K, rho);
    % Translate the parameterizations back
    a = a + perPt(1);
    b = b + perPt(2);
end % end fourierParam