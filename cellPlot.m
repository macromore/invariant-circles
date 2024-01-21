% Script file: cellPlot.m
%
% Purpose: function
% Plot multiple pairs of Fourier parameterizations which are stored in a
% cell array
% 
% Record of revisions:
% Date     Programmer    Description of change
% ====     ==========    =====================
% 09/13/16 David         Original code
% 12/15/16 David         Documenting
%
% Define variables:
% param         - cell array containing pairs of Fourier series
% holdState     - stores the value of plot hold
% 
% Dependencies: 
% biplot.m
%     Fourier.m
% 

function cellPlot(param, color)
    % Check current hold status and turn it on
    holdState = ishold;
    if holdState == 0
        hold on
    end % end if
    % Plot the pairs of Fourier series
    if nargin == 2
        for k = 1:size(param,2)
            biplot(param{1,k},param{2,k}, 2*pi, color);            
        end % end for loop
    else % nargin != 2
        for k = 1:size(param,2)
            biplot(param{1,k},param{2,k}, 2*pi);            
        end % end for loop
    end % end if
    %Turn hold off if it was off to begin with
    if holdState == 0
        hold off
    end % end if
end % end cellPlot