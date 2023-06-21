function phiVector = phiVectorConversion(beta, param)
% convert cell array to vectorized parameterization
%
% Define variables:
% beta      - phase conditions varaible (input)
% param     - parameterization of invariant circles (input)
% phi       - output vector
% K         - period of circles
% 
% Dependencies: 
% N/A
    % put beta in first
    phiVector = beta;
    % determine period of trajectory
    K = size(param,2);
    % add each pair of series
    for k = 1:K
       phiVector = [phiVector; mat(param{1,k}); mat(param{2,k})]; 
    end % end for loop
end % end phiVectorConversion