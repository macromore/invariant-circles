function yval = mapy(x,y,alpha)
% Compute new y-value
% 
% Define variables:
% x      - x coordinate (input)
% y      - y coordinate (input)
% alpha  - Henon map parameter (input)
% yval   - new y coordinate
% 
% Dependencies: 
% N/A
    yval = x*sin(alpha) + y*cos(alpha) - x*x*cos(alpha);
end % end mapy