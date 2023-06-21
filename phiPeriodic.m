%% Script file: phiPeriodicOneRotation.m
%
% Purpose: function
% phi map with one rotation between last circle and first
% 
% Record of revisions:
% Date     Programmer    Description of change
% ============================================
% 09/13/16 David         Original code
% 12/20/16 David         Documenting
% 06/25/20 David         Updated style
%
% Define variables:
% beta      - phase conditions varaible (input)
% param     - parameterization of invariant circles (input)
% alpha     - Henon map parameter (input)
% rho       - rotation number of invariant circle (input)
% phase     - phase condition (input)
% betanext  - output next beta value
% K         - period of circles
% paramnext - output next parametrization
% 
% Dependencies: 
% Fourier.m
% mapx.m
% mapy.m
%% Begin the function
function [betanext, paramnext] = phiPeriodic...
        (beta, param, alpha, rho, phase)
    % copmute next beta value
    betanext = evaluate(param{2,1},0) - phase; % - beta;
    % Compute period
    K = size(param,2);
    % Set up cell array for output
    paramnext = cell(2,K);
    % Apply map to parameterizations first K-1 circles
    for k = 1:K-1
        paramnext{1,k} = mapx(param{1,k}, param{2,k},alpha) - param{1,k+1};
        paramnext{2,k} = mapy(param{1,k}, param{2,k},alpha) - param{2,k+1};
    end % end for loop
    % Apply map to last circle, which includes rotation
    paramnext{1,K} = mapx(param{1,K}, param{2,K},alpha) - ...
        (1+beta)*rotation(param{1,1},rho*K);
    paramnext{2,K} = mapy(param{1,K}, param{2,K},alpha) - ...
        (1+beta)*rotation(param{2,1}, rho*K);
end % end phiPeriodic