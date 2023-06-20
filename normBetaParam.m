function val = normBetaParam(beta, param, nu)
% Compute norm of beta and parmeterization
%
% Define variables:
% beta      - phase condition variable (input)
% param     - parameterization (input)
% nu        - Banach space parameter (input)
% val       - output norm
% 
% Dependencies: 
% norm.m
% 
% Notes: Need to clean up the rotation number computaton to be done first
% and passed to subsequent functions.
    val = zeros(2,size(param,2)+1);
    % val will return max of l_1,nu norms
    val(end) = abs(beta);
    % add up norms of the parameterization pairs
    for k = 1:size(param,2)
       val(1,k) = l1Norm(param{1,k}, nu);
       val(2,k) = l1Norm(param{2,k}, nu); 
    end % end for loop
    val = max(max(val));
end % end normBetaParam