function logCoeffPlot( param, color )
% Plot the log of the coefficents of a fourier series
% 
% Define variables:
% param         - parameterizations of tor(input)
% color         - desired color of plot (input)
% 
% Dependencies: 
% N/A
    if nargin == 1
        color = 'b';
    end % end if
    if isa(param, 'cell')
        maxParam = cell(2, size(param,2));
        for k = 1:2
            for j = 1:size(maxParam,2)
                maxParam{k,j} = matrix(param{k,j});
            end % end for loop
        end % end for loop
        N = length(param{1});
        maxParam = cat(2, maxParam{:});
        maxParam = max(maxParam,[],2);
    else % if it is just a fourier series
        N = length(param);
        maxParam = mat(param);
    end % end if
    plot(-N:N, log(abs(maxParam)), color);
    xlabel('Fourier Coefficeint')
    ylabel('Log_{10} Modulus')
    set(gcf, 'Name', 'Log Plot of Coefficeints')
    axis tight
end % end logCoeffPlot