function traj = trajectoryPlot(P,alpha,numPoints,color)
% Plot trajectory
%
% Define variables:
% P         - point to generate trajectory (input)
% alpha     - area preserving Henon map parameter (input)
% numPoints - length of trajectory to plot (input)
% color     - color of plot (input)
% 
% Dependencies: 
% pointTrajectory.m
    % set up defaults
    if nargin == 1
        alpha = pi/2;
        numPoints = 500;
        color = 'b';    
    elseif nargin == 2
        numPoints = 500;
        color = 'b';
    elseif nargin == 3
        color = 'b';
    end % end if
    % compute trajectory
    orbit = pointTrajectory(P, alpha, numPoints);
    % check hold and turn on
    if size(findall(0,'type','figure')) ~= 0
        ty = ishold;
    else % if there is only one figure
        ty = 0;
    end %end if
    if ty == 0
        figure
        hold on
    end % end if
    % Plot trajectory
    traj = plot(orbit(1, 2:end), orbit(2,2:end), strcat(color,'.'));
    if ty == 0
        hold off
    end % end if
end % end trajectoryPlot