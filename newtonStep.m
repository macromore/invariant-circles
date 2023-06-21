function [betanext, paramnext] = newtonStep(beta, param, alpha, rho, phase)
% Preform a newton step
%
% Define variables:
% beta      - phase conditions varaible (input)
% param     - parameterization of invariant circles (input)
% alpha     - Henon map parameter (input)
% rho       - rotation number of invariant circle (input)
% phase     - phase condition (input)
% x         - vectorized parameterization
% betanext  - output next beta value
% paramnext - output next parametrization
% betatemp  - image of previous beta
% paramtemp - image of previous parameterization
% xnext     - vector version of output
% N         - total length of Fourier series
% 
% Dependencies: 
% phiVectorConversion.m
% phiPeriodic.m
% dPhi.m
% Fourier.m
% phiCellConversion.m
    % convert input to a vector
    x = phiVectorConversion(beta, param);
    % get the phi image of input
    [betatemp,paramtemp] = phiPeriodic(beta, param, alpha, rho, phase);
    % preform a newton step
%     s = toc;
%     if 2*length(param{1})*size(param,2) > 14000
%     if length(param{1})*size(param,2) > 2500
        xnext = x- dPhiSparse(beta, param, alpha, rho)\...
            phiVectorConversion(betatemp,paramtemp);
%         fprintf('Sparse!\n')
%     else
%         xnext = x - dPhi(beta, param, alpha, rho)\...
%             phiVectorConversion(betatemp,paramtemp);
%     end
%         t = toc - s;
%         fprintf('CPU inverse: %g\n',t)
%         s = toc;
%     else
%         xnext = gpuArray(x) - gpuArray(dPhi(beta, param, alpha, rho))...
%             \ gpuArray(phiVectorConversion(betatemp,paramtemp));
%         xnext = gather(xnext);
%         t = toc - s;
%         fprintf('GPU inverse: %g\n',t)
%     end
    % Compute lenght of each series
    N = length(param{1,1});
    % Convert back to cell format
    [betanext, paramnext] = phiCellConversion(xnext, N);
end % end newtonStep