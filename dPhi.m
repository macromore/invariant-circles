%% Script file: dPhi.m
%
% Purpose: function
% Compute the Freche derivative of the conjugacy map 
% 
% Record of revisions:
% Date     Programmer    Description of change
% ============================================
% 09/13/16 David         Original code
% 12/15/16 David         Documenting
% 06/25/20 David         Updated style
%
% Define variables:
% beta      - phase condition (input)
% param     - cell array of parameterizations (input)
% alpha     - parameter for henon map (input)
% rho       - computed rotation number (input)
% N         - number of Fourier modes (derived)
% K         - number of periodic circles (derived)
% n         - length of Fouerier series (derived)
% vfill     - vertical zeros
% hfill     - horizontal zeros
% fill      - square zeros
% dPhiCell  - cell array to hold matrices
% rotMatrix - derivative matrix of a rotaion
% dPhiMatrix- the freche derivative
% 
% Dependencies: 
% Fourier.m
% 
% Notes:
% Should I compute sin(alpha) and cos(alpha) only once?
%% Begin function
function dPhiMatrix = dPhi(beta, param, alpha, rho)
    % Initialize variables (most of them)
    N = length(param{1,1});
    K = size(param,2);
    n = 2*N+1;
    vfill = zeros(n,1);
    hfill = zeros(1,n);
    fill = zeros(n);
    dPhiCell = cell(2*K+1, 2*K+1);
    dPhiCell{1,1} = 0;
    % Fill the square protions of the derivative with zeros
    for k = 2:2*K+1
        for j = 2:2*K+1
            dPhiCell{k,j} = fill;
        end % end for loop
    end % end for loop
    % Fill the horizontal and vertical portions of the derivative with zeros
    for k = 2:2*K+1
        dPhiCell{1,k} = hfill;
        dPhiCell{k,1} = vfill;
    end % end for loop
    % Apply the derivative with respect to beta?
    dPhiCell{1,3} = dPhiCell{1,3} - ones(1,n);
    dPhiCell{2*K,1} = dPhiCell{2*K,1} - mat(rotation(param{1,1},K*rho));
    dPhiCell{2*K+1,1} = dPhiCell{2*K+1,1} - mat(rotation(param{2,1},K*rho));
    % Apply with respect to "next" periodic circle 
    for k = 1:K-1
        dPhiCell{2*k,2*k+2} = dPhiCell{2*k,2*k+2} - eye(n);
        dPhiCell{2*k+1,2*k+3} = dPhiCell{2*k+1,2*k+3} - eye(n);
    end % end for loop
    % Set up rotation matrix
    rotMatrix = zeros(n);
    for k = -N:N
        rotMatrix(N+k+1,N+k+1) = exp(2*pi*1i*(K*rho)*k);
    end % end for loop
    % Apply derivative to final periodic circle
    dPhiCell{2*K,2} = dPhiCell{2*K,2} - (1+beta)*rotMatrix;
    dPhiCell{2*K+1,3} = dPhiCell{2*K+1,3} - (1+beta)*rotMatrix;
    % Fill in derivative for the henon portions
    for k = 1:K
        dPhiCell{2*k,2*k} = dPhiCell{2*k,2*k} + mat(cos(alpha)*eye(n) + ...
            FourierOperator(2*sin(alpha)*param{1,k}));
        dPhiCell{2*k+1,2*k} = dPhiCell{2*k+1,2*k} + mat(sin(alpha)*eye(n) - ...
            FourierOperator(2*cos(alpha)*param{1,k}));
        dPhiCell{2*k,2*k+1} = dPhiCell{2*k,2*k+1} - sin(alpha)*eye(n);
        dPhiCell{2*k+1,2*k+1} = dPhiCell{2*k+1,2*k+1} + cos(alpha)*eye(n);
    end % end for loop
    %Convert the cell array to a matrix
    dPhiMatrix = cell2mat(dPhiCell);
end % end function