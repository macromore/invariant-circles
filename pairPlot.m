function pairPlot(a, b, color)
% Graph a pair of Fourier series
% 
% Define variables:
% a           - x coordinate fourier series (input)
% b           - y coordinate fourier series (input)
% omega       - fourier series frequency (input)
% tArray      - time array for evalutation
% x           - x values at times
% y           - y values at times
% 
% Dependencies: 
% Fourier.m
    if nargin == 2
        color = 0;
    end %end if
    % Set up the time steps
    tArray = 0:1/10000:1;
    % Evaluate the Fourier series at the times
    x = real(evaluate(a,tArray));
    y = real(evaluate(b,tArray));
    % Plot the resulting points
    if color == 0
        plot(x, y);
    else % color != 0
        plot(x,y,color);
    end % end if
end % end function