function val = normPhi(beta, param, alpha, rho, phase, nu)
% apply phi and take the norm of a parameteriation
% 
% Define variables:
% beta      - phase condition variable (input)
% param     - parameterization (input)
% alpha     - Henon map parameter (input)
% rho       - rotation number(input)
% phase     - phase condition (input)
% nu        - Banach space parameter (input)
% newbeta   - newton step computed beta
% newparam  - newton step computed parameterization
% val       - defect under the phi map
% 
% Dependencies: 
% phiPeriodic.m
% normBetaParam.m
    % Apply phi to the parameterization
    [newbeta, newparam] = phiPeriodic(beta, param, alpha, rho, phase);
    % Take norm of image
    val = normBetaParam(newbeta, newparam, nu);
end % end normPhi