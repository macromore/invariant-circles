function sobolevNorm = sobolevNorm(param, sobolevParam)
% Compute the sobolev norm of a set of parameterizations
%
% Define variables:
% param       - cell array of Fourier parameterizations (input)
% s           - H^s (input)
% sobolevNorm - sobolev norm array (output)
% x           - temp variable for matrices
% K           - number of Fourier modes
% sum         - temp variable to hold norm
% c           - holds both values for the parameterizations
% 
% Dependencies: 
% Fourier.m
    if nargin == 1
        sobolevParam = 1;
    end % end if
    sobolevNorm = zeros(size(param,2),1);
    for j = 1:size(param,2)
        x = [mat(param{1,j})';mat(param{2,j})'];
        K = length(param{1,j});
        sum = 0;
        for k = -K:K
            c = x(:,k+K+1);
            sum = sum + (1+abs(k)^2)^sobolevParam*max(abs(c))^2;
        end % end for loop
        sobolevNorm(j) = sqrt(sum);
    end % end for loop
end % end sobolevNorm

