function ppoints = periodicPointFinder(orbits, K, alpha)
% Preform newton's method on an orbit to find periodic point in the center
% for use parameterizting. 
% 
% Define variables:
% orbits  - trajectory to be analized (input)
% K       - Number of periodic points (input)
% alpha   - parameter for henon (input)
% ppoints - periodic points (output)
% rho     - dummy rotation number = 1 
% param   - intermediate variable for periodic points
% beta    - newton parameter
% derr    - difference in error
% err     - netwon error
% err1    - holds new error value
% j       - iteration counter
% 
% Dependencies: 
% Fourier.m
% newtonStep.m
% 
% Notes:
% Put dummy values in for the reqired things
% We use the mean because the boundary does not always converge to the
% correct fixed points
    ppoints = zeros(2,K);
    rho = 1;
    param = cell(2,K);
    beta = 0;
    derr = 1;
    err = 1;
    iterationCount = 1;
    for k = 1:K
       % Use average of orbit as initial guess
       ppoints(1,k) = mean(orbits{k}(1,:));
       ppoints(2,k) = mean(orbits{k}(2,:));
       % Convert to Fourier series
       param{1,k} = Fourier(ppoints(1,k));
       param{2,k} = Fourier(ppoints(2,k));
    end % end for loop
    plot(ppoints(1,:),ppoints(2,:),'.m')
    % Run a newton step with conjugacy check
    if K > 1
        while derr > 10^-10 && iterationCount < 10
             [beta, param] = newtonStep(beta, param, alpha, rho, ...
                 evaluate(param{2,1},0));
             temp = [evaluate(param{1,1},0); evaluate(param{2,1},0)];
             for k = 1:K
                 temp = henon(temp, alpha);
             end % end for loop
             err1 = norm([evaluate(param{1,1},0); ...
                 evaluate(param{2,1},0)] - temp);
             derr = err - err1;
             err = err1;
             iterationCount = iterationCount + 1;
        end % end while loop
        fprintf('Periodic point conjugacy error is %1.15f\n', err)
        for k = 1:K   
           ppoints(1,k) = evaluate(param{1,k},0);
           ppoints(2,k) = evaluate(param{2,k},0);   
        end %end for loop
        %plot(ppoints(1,:),ppoints(2,:),'*')
    end % end if
end % end periodicPointFinder