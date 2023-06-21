% Analytic continuation with Newton and rho finding %%%%%%%%%%%%%%%%%%%%%%%
%clear, close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
% alpha = acos(.24);
% % Period 1
%     P = [0.4,0];
% % Period 5
%    P = [0.5, 0];
%    P = [0.6, 0];
%    P = [0, 0.46];
%    P = [.4, -.4];

% alpha = 3*pi/2;
% % Period 1
%    P = [0.3,0];
% % Period 21
%    P = [-0.676, 0]; % does not appear to work with NC
% % Period 29
%    P = [-.5, -.15];
% % Period 25
%    P = [-.5, -.2];


% alpha = acos(-0.95);
% % Period 1
%     P = [.2,0];
%     P = [0,-.4];
%     P = [-.3, -2]; % Fails, large starting norm
%     P = [.3,-1.8];
%     P = [-.2, -2.5]; % Works large starting norm
    
% % Period 7
%    P = [0,-1.9];
% 
% % Period 12
%     % P = [0,-2.67];
% 
% % Period 120
%     % P = [0; -2.65];
% Period 400 % Not real
%     % P = [-.6323, .2671];
% 

% % Figure 1.38
% alpha = acos(.8);
% % Period 1
% P = [0.05, 0];
% % Period 26 - doesn't work maybe some tweaking...
% P = [0.52835,0];
% % Period 26*like 4
% P = [0.528355,0];

% % Figure 1.39
% alpha = acos(.4);
% % % Period 1 : Center
% % P = [.5,0];
% % % Period 6
% % P = [.6, .2];
% % % Period 1 : Again
% % P = [.7, 0];
% % % Period 25
% % P = [.73, .43];
% % Period 1 - Large
% % P = [.7, .455];
% % Period 19
% % P = [.74, .48];
% % Period 1 - doesn't grow
% % P = [.7,.59];
% % Test
% P = [.74,0];

% Figure 1.40
alpha = acos(.24);
% % Period 1
% P = [.4,0];
% % Period 5
% P = [.55,0];
% % Period 5 (didn't stop)
% P = [.6, -.06];
% Test
P = [.75,0];

% Set up variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Banach Space Parameter
nu = 1.1;
% Initial number of modes, and how many to add
initialmodes = 50;
modestep = 50;
modes = initialmodes;
maxModes = 400;
% Number of points for plotting and computing initial parameterization
numpoints = 10000;
% Set total time to zero
t = 0;
% Largest value allowed to change rho
maxRhoStep = 10^-4;
% Do we need to replace rho with 1-rho
rhoFlip = 1;
% Do we need to subtract or add to rho
delRhoFlip = 1; % Sign flip on delRho
% How small should be go on the steps
maxPower = 5;
% Error limit on rho
rhoLimit = 14;
% Error limit in Newton
errorLimit = 14;
% Error limit in truncation
tailLimit = 16;
% Sobolev space H^m max
mMax = 10;
% Attempt limit
maxAttempts = 200;

% Command line formatting and figure hold %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long
hold on
attemptNumber = 1;
successNumber = 0;

% Ask for a file name for the diary output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = input('Name folder for output? (can be root/folder, etc.) ','s');
mkdir(fname);
diary(strcat('./',fname,'/output.txt'))
successFile = fopen(strcat('./', fname, '/successes.txt'), 'w');

% Header information for the file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Parameters: \n')
fprintf('P =                      [ %f, %f ]\n', P(1),P(2))
fprintf('alpha =                  %.15g                      probably acos(something)\n', alpha)
fprintf('nu =                     %.15g\n', nu)
fprintf('sobolev max space m =    %d\n', mMax)
fprintf('initial number of modes  %d\n', initialmodes)
fprintf('Mode step used in newton %d\n', modestep)
fprintf('number of points         %d\n', numpoints)
fprintf('maximum rho step         %.15g\n', maxRhoStep)
fprintf('minimum rho step         %.15g\n', maxRhoStep*10^-maxPower)
fprintf('max tail size            %.15g\n', 10^-tailLimit)
fprintf('max error                %.15g\n', 10^-errorLimit)
fprintf('rho difference limit     %.15g \n', 10^-rhoLimit)
fprintf('rho flip                 %d\n', rhoFlip)
fprintf('direction of rho change  %d\n', delRhoFlip)
fprintf(successFile, '\n'); 
fprintf(successFile, 'Parameters: \n');
fprintf(successFile, 'P =                      [ %f, %f ]\n', P(1),P(2));
fprintf(successFile, 'alpha =                  %.15g                      probably acos(something)\n', alpha);
fprintf(successFile, 'nu =                     %.15g\n', nu);
fprintf(successFile, 'sobolev space max m =    %d\n', mMax);
fprintf(successFile, 'initial number of modes  %d\n', initialmodes);
fprintf(successFile, 'Mode step used in newton %d\n', modestep);
fprintf(successFile, 'number of points         %d\n', numpoints);
fprintf(successFile, 'maximum rho step         %.15g\n', maxRhoStep);
fprintf(successFile, 'minimum rho step         %.15g\n', maxRhoStep*10^-maxPower);
fprintf(successFile, 'max tail size            %.15g\n', 10^-tailLimit);
fprintf(successFile, 'max error                %.15g\n', 10^-errorLimit);
fprintf(successFile, 'rho difference limit     %.15g \n', 10^-rhoLimit);
fprintf(successFile, 'rho flip                 %d\n', rhoFlip);
fprintf(successFile, 'direction of rho change  %d\n', delRhoFlip);
fprintf(successFile, '\n');

% Compute and plot trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Computing the trajectory... \n    ')
tic
trajectory = pointTrajectory(P, alpha, numpoints);
toc
t = t + toc;
fprintf('\n')

% Plot trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(trajectory(1,:),trajectory(2,:),'.r');
%savefig(strcat('./',fname,'/fig1_trajectory'));

% Ask the number of tori and separate the orbits %%%%%%%%%%%%%%%%%%%%%%%%%%
K = input('How many periodic tori? ');
fprintf('\nSeparate the trajecotory... \n    '    )
tic
tori = trajectorySeparator(trajectory, K);
toc
t = t + toc;
fprintf('\n')

% Find periodic points of Tori %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Finding periodic points ... \n    ')
tic
perPts = periodicPointFinder(tori, K, alpha);
toc
t = t + toc;
fprintf('\n')

% Plot periodic points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(perPts(1,:),perPts(2,:),'.k')
%savefig(strcat('./',fname,'/fig2_perPts'));

% Find rotation number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Compute the rotation number... \n')
tic
i = 1;
diffRho = 1;

while diffRho > 10^-rhoLimit && numpoints < 10^8
    numpoints = K*100*2^i; % Start with 100 points per tori and double it each iteration
    fprintf('    Number of points per tori         %d\n    ', 100*2^i)
    tic
    trajectory2 = pointTrajectory(P, alpha, numpoints);
    tori2 = trajectorySeparator(trajectory2, K);
    toc
    rho(i) = weightedBirkoffRotationNumber(tori2, K, perPts);
    fprintf('    Rho approximates as %.15g. \n    ',rho(i))
    toc
    if i > 1
        diffRho = abs(rho(i) - rho(i-1));
        fprintf('    Difference between last rho and this rho: %.15d\n', abs(rho(i)-rho(i-1)))
        fprintf('    Machine epsilon is %.15d. \n\n', eps(rho(i)))
    else
        fprintf('\n')
    end
    i = i + 1;
end
if numpoints >= 10^8
    fprintf('Rho accuracy questionable.\n')
end
fprintf(successFile, 'The rotation number is approximately %.15g.\n    ', rho(end));
fprintf('The rotation number is approximately %.15g.\n    ', rho(end))
rho = rho(end);
toc
t = t + toc;
fprintf('\n')

% Initialize parameterization and pad with zeros %%%%%%%%%%%%%%%%%%%%%%%%%%
param = cell(2,K);
fprintf('Compute the initial Fourier modes with %d initial modes... \n    ', initialmodes)
tic
for i = 1:K
   [param{1,i}, param{2,i}] = fourierParam(tori{i}, initialmodes, K, rho, perPts(:,i)); 
   param{1,i} = truncate(param{1,i},modes);
   param{2,i} = truncate(param{2,i},modes);
end
toc
t = t + toc;
fprintf('\n')

% Plot initial parameterization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellPeriodicPlot(param, 'm');
% savefig(strcat('./',fname,'/fig_3initial_param'));

% rhoFlip = input('Do we need 1-rho for the newton-like operator (0 or 1)? ');
if rhoFlip == 1 % rhoFlip(1) == '1'
    rho = 1 - rho;
    fprintf('    Flipped rho is %.15g.\n\n', rho)
else
    fprintf('\n')
end

% Measure the sequence space error of the initial parameterization %%%%%%%%
fprintf('Compute initial parameterization error with nu = %1.1f ... \n    ', nu)
tic
beta = 10^-K;
phase = evaluate(param{2,1},0,2*pi);
error = normPhi(beta, param, alpha, rho, phase, nu);
fprintf('Error = %.15g \n     ', error)
toc
t = t + toc;
fprintf('\n')

% % % Set up tolerance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % J = input('What precision to stop at, 10^-J? (suggest about 10^-14) ');
% % fprintf('\n')

iteration = 0;
tic

% Initial Newton-like operation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Carry out newton-like method...\n')
while error > 10^-errorLimit && iteration < 10
    [beta, param] = newtonStep(beta, param, alpha, rho, phase);
    iteration = iteration + 1;
    error = normPhi(beta, param, alpha, rho, phase, nu);
    fprintf('    Iteration %d, error = %.15g \n     ', iteration, error)
    toc
end
fprintf('\n')

% Sobolev norms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sobolevGrid = zeros(mMax,1);
    
% Compute initial sobolev norm and tail value %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ev = endVals(param{1});
for m = 1:mMax
    maxSobolevNorm = max(sobolevNorm(param,m));
    fprintf(successFile, 'max Sobolev m = %d Norm is %.15g \n', m, maxSobolevNorm);
    fprintf('max Sobolev m = %d Norm is %.15g \n', m, maxSobolevNorm)
    sobolevGrid(m, 1) = maxSobolevNorm;
end
fprintf(successFile, 'last coefficient absolute value is %.15g \n', ev);
fprintf('last coefficient absolute value is %.15g \n', ev)
t = t + toc;
fprintf('Total time elapsed %f seconds. \n', t)

% Plot initial newton output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellPeriodicPlot(param, 'k')
legend('Trajectory','Initial Guess','Approximation')
savefig(strcat('./',fname,'/fig4_first_param'));

% Plot initial newton output on new figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
cellPeriodicPlot(param, 'k')

% Plot and save initial coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
logCoeffPlot(param{1})
savefig(strcat('./',fname,'/iniitialCoeffs'));
% Analytic continuation time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_original = rho;
param_original = param;

% Confirem contuation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input('\nContinue? (Ctrl-C to stop): ', 's');
i = 0;
reRun = 0;

% Main loop, shrinking rho size on unsuccessful attempts %%%%%%%%%%%%%%%%%%
while i <= maxPower && attemptNumber < maxAttempts
    del_rho = ((-1)^delRhoFlip)*maxRhoStep*2^-i;
    fprintf('\nUsing rho step of %.15g \n', del_rho)
    error = 0;

    % Iteration loop, attempting constant rho step size %%%%%%%%%%%%%%%%%%%
    while error < 10^-errorLimit && attemptNumber < maxAttempts
        
        % Make rho smaller %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rho = rho - del_rho;
        if reRun == 0
            fprintf('New rho =                               %.15f\n\n', rho)
        else
            fprintf('Rho =                                   %.15f\n\n', rho)
        end
        
        % Set up newton-like method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        iteration = 0;
        if reRun == 0
            param_old = param;
        end
        fprintf('Start the new newton-like method using %d modes and nu = %1.1f ... \n', modes, nu)
        tic
        beta = 10^-K;
        phase = evaluate(param{2,1},0,2*pi);
        error = normPhi(beta, param, alpha, rho, phase, nu);
        fprintf('    Iteration 0, error = %.15g \n     ', error)
        toc
        t = t + toc;

        % Carry out newton-like operator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        while error > 10^-errorLimit && error < 50 && iteration < 8
            [beta, param] = newtonStep(beta, param, alpha, rho, phase);
            iteration = iteration + 1;
            error = normPhi(beta, param, alpha, rho, phase, nu);
            fprintf('    Iteration %d, error = %.15g \n     ', iteration, error)
            toc
        end
        fprintf('\n')        
        % Compute output metrics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ev = max(max(cellfun(@endVals, param)));
        fprintf('The end coefficients absolute value is: %.15g \n', ev)
        t = t + toc;
        fprintf('    Total time elapsed %f seconds. \n', t)
        % Plot coefficient absolute values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(3)
        logCoeffPlot(param{1})
        
        % Decide next move: continue, rerun with more modes, %%%%%%%%%%%%%%
        % rerun with smaller rhostep, stop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Accept new parameterization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ev < 10^-tailLimit && error < 10^-errorLimit
            figure(2)
            cellPeriodicPlot(param)
            
            % Record success            
            successNumber = successNumber + 1;
            fprintf(successFile, 'Success number %d \n', successNumber);
            fprintf(successFile, 'Attempt number %d, successful!\n', attemptNumber);
            fprintf(successFile, '    number of modes    %d\n', modes);
            fprintf(successFile, '    error              %.15g\n', error);
            fprintf(successFile, '    tail value         %.15g\n', ev);
            fprintf(successFile, '    rho                %f\n', rho);
            
            % Compute Sobolev norms
            tempSobVect = zeros(mMax, 1);
            for m = 1:mMax
                maxSobolevNorm = max(sobolevNorm(param,m));
                fprintf('The max Sobolev m = %d Norm is:                %.15g \n', m, maxSobolevNorm)
                fprintf(successFile, '    sobolev norm m = %d   %.15g\n', m, maxSobolevNorm);
                tempSobVect(m, 1) = maxSobolevNorm;
            end
            sobolevGrid = [sobolevGrid, tempSobVect];
            fprintf(successFile, '\n');
            fprintf('Attempt number %d, successful!\n', attemptNumber)
            fprintf('Success number %d!\n', successNumber)
            attemptNumber = attemptNumber + 1;
            
            % Truncate back to original number of modes and exit loop
            modes = initialmodes;
            for k = 1:K
                param{1,k} = truncate(param{1,k},modes);
                param{2,k} = truncate(param{2,k},modes);
            end
            i = 0;
            error = 1;
            reRun = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Rerun with more modes or move on if too many modes %%%%%%%%%%%%%%
        elseif ev > 10^-tailLimit
            rho = rho + del_rho;
            modes = modes + modestep;
            for k = 1:K
                param{1,k} = truncate(param{1,k},modes);
                param{2,k} = truncate(param{2,k},modes);
            end
            if modes <= maxModes
                error = 0;
                fprintf('\nMore modes required, rerunning Newton with %d modes.\n', modes)
                reRun = 1;
            else
                fprintf('\nMax modes exceeded, trying again with smaller rho step.\n')
                modes = initialmodes;
                param = param_old;
                for k = 1:K
                    param{1,k} = truncate(param{1,k},modes);
                    param{2,k} = truncate(param{2,k},modes);
                end
                reRun = 0;
                rho = rho + del_rho;
                i = i + 1;
                fprintf('\nAttempt number %d, failed!\n', attemptNumber)
                attemptNumber = attemptNumber + 1;
                error = errorLimit;
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Try smaller stepsize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            modes = initialmodes;
            param = param_old;
            for k = 1:K
                param{1,k} = truncate(param{1,k},modes);
                param{2,k} = truncate(param{2,k},modes);
            end
            reRun = 0;
            rho = rho + del_rho;
            i = i + 1;
            fprintf('\nAttempt number %d, failed!\n', attemptNumber)
            attemptNumber = attemptNumber + 1;
        end
    end
end

% Plot the sobolev norms
figure(4)
if size(sobolevGrid,2) == 1
    sobolevGrid = [sobolevGrid, zeros(mMax,1)];
end
mesh(sobolevGrid)

% Save figures and close stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
savefig(strcat('./',fname,'/sobolevNorms'))

figure(3)
savefig(strcat('./',fname,'/finalCoeffs'))

figure(2)
savefig(strcat('./',fname,'/fig4_ac_final_param'))
fclose(successFile);
diary off