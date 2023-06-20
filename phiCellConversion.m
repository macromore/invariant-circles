function [beta, param] = phiCellConversion(vect, N)
% convert vectorized parameter to a cell array
%
% Define variables:
% vect      - vectorized parameteriztion (input)
% N         - total length of fourier series (input)
% beta      - phase condition variable (input)
% vec       - truncated and reshaped vect
% K         - number of periodic circles
% param     - cell version of vect
% 
% Dependencies: 
% Fourier.m
    % Take first element of vector as beta
    beta = vect(1);
    % Slice off first element so only Fourier series remain
    vec = vect(2:end);
    % Compute number of tori
    K = length(vec)/(2*(2*N+1));
    % Rearrange vec into rows of series
    vec = reshape(vec,2*N+1,2*K);
    % Set up parameter cell arrray
    param = cell(2,K);
    % Fill in cell array
    for k = 1:K
        param{1,k} = Fourier(vec(:,2*k-1));
        param{2,k} = Fourier(vec(:,2*k));
    end % end for loop
end % end phiCellConversion