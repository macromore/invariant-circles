function xval = mapx(x,y,alpha)
% Compute new x-value
%
% Define variables:
% x      - x coordinate (input)
% y      - y coordinate (input)
% alpha  - Henon map parameter (input)
% xval   - new x coordinate
% 
% Dependencies: 
% N/A
    xval = x * cos(alpha) - y * sin(alpha) + x * x * sin(alpha);
end % end mapx