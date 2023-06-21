function tori = trajectorySeparator(trajectory, K)
% Separate out a trajectory assuming they are fixed tori of the K-fold
% composition.
%
% Define variables:
% trajectory - the points on an orbit (input)
% K          - period of suspected invariant tori (input)
% tori       - separated trajector (output)
% 
% Dependencies: 
% N/A
    tori = cell(1,K);
    for j = 1:K
       tori{1,j} = trajectory(:,j:K:end);
    end % end for loop
end % end trajectorySeparator