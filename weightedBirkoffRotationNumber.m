%% Script file: weightedBirkhoffRotationNumber.m
% 
% Purpose: function
% Given a periodic trajectory, compute a rotation number using weighted
% averages
% 
% Record of revisions:
% Date     Programmer    Description of change
% ============================================
% 12/15/16 David         Original code
% 12/19/16 David         Documenting
% 12/22/16 David         Changed translation to use first iterate
% 05/26/17 David         Set up to accept periodic point instead of average
% 06/25/20 David         Updated style
%
% Define variables:
% orbit         - initial data, the trajectory (input)
% K             - number of closed orbits (input)
% rho           - computed rotation number
% avex          - average of the x trajectory
% avey          - average of the y trajectory
% numPoints     - number of iterates to use
% args          - angle coordinates of points
% argDiff       - change in angle between two points
% totalWeight   - denominator for weight function
% weightValue   - numerator for weight function
% partialSum    - partial sum of weighted average
% 
% Dependencies: 
% N/A
%% Begin the Function
function rho = weightedBirkoffRotationNumber(orbits, K, perPts)
    % set rho to zero as a trigger for recursion
    rho = 0;
    if K ~= 1
       rho = zeros(1,K);
       for k = 1:K
          rho(k) = weightedBirkoffRotationNumber(orbits{k}, 1, perPts(:,k));
       end % end for loop
       rho = median(rho)/K;
    end % end if
    % if we haven't computed a rotation number yet, compute rotation number
    if rho == 0
        if isa(orbits,'cell')
           orbits = orbits{1}; 
        end % end if
        % translate torus to origin
        orbits(1,:) = orbits(1,:) - perPts(1);
        orbits(2,:) = orbits(2,:) - perPts(2);
        args = atan2(orbits(2,:),orbits(1,:))/(2*pi);
        argDiff = mod(diff(args),1);
        % Compute the weighted sum
        weights = vectorWeight(length(orbits));
        rho = sum(weights(1,2:end) .* argDiff);
    end % end if
end % end weighted BirkhoffRotationNumber