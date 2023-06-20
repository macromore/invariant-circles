function phasespacePlot(alpha,numTrajectories, numPoints, color)
% Plot some trajectories
%
% Define variables:
% alpha           - henon parameter (input)
% numTrajectories - number of trajectories to plot (input)
% numPoints       - number of points per trajectory (input)
% color           - color of the plot (input)
% P               - linear space for the plot
% scale           - size of the plot
% ty              - checks the hold state
% 
% Dependencies: 
% N/A
     scale = 4;
     if nargin == 0
        alpha = pi/2; 
        numTrajectories = 100;
        numPoints = 200;
        color = 'k';
    elseif nargin == 1
        numTrajectories = 100;
        numPoints = 100;
        color = 'k';
    elseif nargin == 2
        numPoints = 100;
        color = 'k';
    elseif nargin == 3
        color = 'k';
    end % end if
    P = linspace(-1.4*scale, 1.4*scale, numTrajectories);
    scale = 1;
    if size(findall(0,'type','figure')) ~= 0
        ty = ishold;
    else % if there is only one figure
        ty = 0;
    end %end if
    if ty == 0
        figure
        hold on
    end % end if
    for n = 1:numTrajectories
        trajectoryPlot([P(n),P(n)], alpha, numPoints, color);
    end % end for loop
    axis([-scale, scale, -scale, scale])
    if ty == 0
        hold off
    end % end if
end % end phasespacePlot