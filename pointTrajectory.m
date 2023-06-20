function trajectory = pointTrajectory(P, alpha, numPoints)
% Compute a point trajectory
%
% Define variables:
% P          - starting point of trajectory (input)
% alpha      - henon parameter (input)
% numPoints  - number of points in the trajectory (input)
% trajectory - computed trajectory (output)
% cosa       - cos(alpha)
% sina       - sin(alpha)
% 
% Dependencies: 
% N/A
    trajectory = zeros(2,numPoints);
    trajectory(:,1) = P;
    cosa = cos(alpha);
    sina = sin(alpha);
    for k = 2:numPoints
        trajectory(:,k) = [trajectory(1,k-1)*cosa - trajectory(2,k-1)*...
            sina + trajectory(1,k-1)*trajectory(1,k-1)*sina; ...
            trajectory(1,k-1)*sina + trajectory(2,k-1)*cosa - ...
            trajectory(1,k-1)*trajectory(1,k-1)*cosa];
    end % end for loop
end % end pointTrajectory