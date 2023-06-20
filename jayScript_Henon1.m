% Preform analytic continuation and save the pertinent information.
% 
% Notes: 
% This is a fairly complicated piece of script (to me), so there has been
% an attempt to document every step of the procedure. 
% Please note, only the "Map and Trajectory Parameters" and "Set up
% parameter variables" sections need to be modified to run a typical
% computation.
%% Clean MatLab %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearAll
%% Map and Trajectory Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other ideas for starting points are in parameters.m
% alpha = acos(-0.95);
%%%% Period 1 % Use 10^-3 % To Grow delRhoFlip = 1
% initialP = [0,.1]; 
%%%% Period 7 % Use 10^-4 % To Grow delRhoFlip = 0
% initialP = [.3,-1.8];
% initialP = [1.5,2];
%%%% Period 1 Again % Use 10^-4? % To grow delRhoFlip = 1
% initialP = [1.9,1.6]; 
% initialP = [1.8, 1.6];
% initialP = [2.1, 2.4];
% Pretty clearly something between these two
%%%% Period 19 %  To grow delRhoFlip = 0
% initialP = [2.1,3.7];
%%%% Period 1
% initialP = [2.2,2.7]; % To grow delRhoFlip = 1
%%%% Period 12 % Use 10^-4 % To grow delRhoFlip = 0
% initialP = [0,-2.7];
%%%% Period 120 % Use 10^-6? 
% initialP = [0; -2.65];
%% Set up parameter variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Banach Space Parameter
nu = 1.1;
% Initial number of modes, how many to add, and maximum number of modes
initialModes = 10;       %Compute with Yorke/Sander
modeStep = 50;           %If Newton fails
maxModes = 310;          %Max Number of Modes
% Number of points for plotting and computing initial parameterization
numPoints = 1e4;
% Largest value allowed to change rho
maxRhoStep = 10^-3;
% Do we need to replace rho with 1-rho
rhoFlip = 1;
% Do we need to subtract or add to rho (-1^delRhoFlip)
delRhoFlip = 1; 
% How small should be go on the steps 2^-maxpower
maxPower = 2;
% Error limit on rho
rhoLimit = 16;       %10^(-rhoLimit)
% Error limit in Newton
errorLimit = 14;     %10^(-errorLimit)
% Error limit in tail of truncated series
tailLimit = 14;
% Sobolev space H^m max
sobolevMax = 10;
% Conjugacy check max
conjMax = 1000;
% Attempt limit
maxAttempts = 5;
%% Example 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = acos(0.24);
initialP = [0.4,0];
initialModes = 30;
numPoints = 1e2;
maxAttempts = 0;
%% Example 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alpha = acos(0.24);
% initialP = [0.5,0];
% initialModes = 30;
% numPoints = 1e3;
% maxAttempts = 0;
%% Example 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alpha = acos(-0.95);
% initialP = [0,-2.65];
% initialModes = 15;
% numPoints = 12000;
% maxAttempts = 0;
%% Plot color %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color = 'g';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% For general use, nothing below here needs to be changed.                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up Loop Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
attemptNumber = 1;           % Counter for attempts
successNumber = 0;           % Counter for successes
currentModes = initialModes; % numer of modes presently being used
%% Ask for a file name for the diary output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = input('Name folder for output? (can be root/folder, etc.) ','s');
mkdir(fname);
diary(strcat('./',fname,'/output.txt'))
successFile = fopen(strcat('./', fname, '/successes.txt'), 'w');
%% Header information for the file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Parameters: \n')
fprintf(['P =', repmat(' ',1,22),'[ %f, %f ]\n'], initialP(1),initialP(2))
fprintf(['alpha =', repmat(' ',1,18), '%.15g', repmat(' ',1,10),...
    'probably acos(something)\n'], alpha)
fprintf('nu =                     %.15g\n', nu)
fprintf('sobolev max space m =    %d\n', sobolevMax)
fprintf('initial number of modes  %d\n', initialModes)
fprintf('Mode step                %d\n', modeStep)
fprintf('points per circle        %d\n', numPoints)
fprintf('maximum rho step         %.15g\n', maxRhoStep)
fprintf('minimum rho step         %.15g\n', maxRhoStep*2^-maxPower)
fprintf('max tail size            %.15g\n', 10^-tailLimit)
fprintf('max error                %.15g\n', 10^-errorLimit)
fprintf('Conjugacy iterations     %d\n', conjMax)
fprintf('rho difference limit     %.15g \n', 10^-rhoLimit)
fprintf('rho flip                 %d\n', rhoFlip)
fprintf('direction of rho change  %d\n\n', delRhoFlip)
fprintf(successFile, '\n'); 
fprintf(successFile, 'Parameters: \n');
fprintf(successFile, ['P =', repmat(' ',1,22),'[ %f, %f ]\n'], ...
    initialP(1),initialP(2));
fprintf(successFile, ['alpha =', repmat(' ',1,18), '%.15g', ...
    repmat(' ',1,20), 'probably acos(something)\n'], alpha);
fprintf(successFile, 'nu =                     %.15g\n', nu);
fprintf(successFile, 'sobolev space max m =    %d\n', sobolevMax);
fprintf(successFile, 'initial number of modes  %d\n', initialModes);
fprintf(successFile, 'Mode step                %d\n', modeStep);
fprintf(successFile, 'points per circle        %d\n', numPoints);
fprintf(successFile, 'maximum rho step         %.15g\n', maxRhoStep);
fprintf(successFile, 'minimum rho step         %.15g\n', ...
    maxRhoStep*2^-maxPower);
fprintf(successFile, 'max tail size            %.15g\n', 10^-tailLimit);
fprintf(successFile, 'max error                %.15g\n', 10^-errorLimit);
fprintf(successFile, 'Conjugacy iterations     %d\n', conjMax);
fprintf(successFile, 'rho difference limit     %.15g \n', 10^-rhoLimit);
fprintf(successFile, 'rho flip                 %d\n', rhoFlip);
fprintf(successFile, 'direction of rho change  %d\n', delRhoFlip);
fprintf(successFile, '\n');
%% Compute trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Computing the trajectory... \n  ')
tic
trajectory = pointTrajectory(initialP, alpha, numPoints);
toc
fprintf('\n')
%% Plot trajectory and figure hold %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figTrajectory = figure(1);
set(figTrajectory, 'Units', 'Normalized', 'OuterPosition', ...
    [.05 .5 .25 .4]);
set(figTrajectory, 'Name', 'Trajectory and Initial Parameterization');
hold on
phasespacePlot(alpha)
plot(trajectory(1,:),trajectory(2,:),'.g');
axis([-4 4 -4 4]) % Do I need to do this?

times = linspace(0, 99, 100)
figure 
hold on
plot(times, trajectory(1, :), 'b')
plot(times, trajectory(2, :), 'g')


%% Ask the number of tori and separate the orbits %%%%%%%%%%%%%%%%%%%%%%%%%
K = 1 %input('How many periodic tori? ');
fprintf('\nSeparate the trajecotory... \n  ')
tori = trajectorySeparator(trajectory, K);
toc
fprintf('\n')
%% Find periodic points of Tori %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Finding periodic points ... \n  ')
perPts = periodicPointFinder(tori, K, alpha);
toc
fprintf('\n')
plot(perPts(1,:),perPts(2,:),'.k')
%% Find rotation number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Compute the rotation number... \n')
rhoMultiplier = 1;
diffRho = 1;
numPointsRho = 0;
while diffRho > 2.224*10^-rhoLimit && numPointsRho < 10^10
    % Start with 1000 points per tori and add 1000 each iteration
    numPointsRho = K*100*rhoMultiplier; 
    fprintf('Number of points per tori         %d\n', 1000*rhoMultiplier)
    trajectory2 = pointTrajectory(initialP, alpha, numPointsRho);
    tori2 = trajectorySeparator(trajectory2, K);
    rho(rhoMultiplier) = weightedBirkoffRotationNumber(tori2, K, perPts);
    fprintf('  Rho approximates as %.15g. \n',rho(rhoMultiplier))
    if rhoMultiplier > 1
        diffRho = abs(rho(rhoMultiplier) - rho(rhoMultiplier-1))
        fprintf('  Delta rho:          %e \n'...
            , abs(rho(rhoMultiplier)-rho(rhoMultiplier-1)))
    end % end if
    fprintf('  Machine epsilon is  %e \n  ', eps(rho(rhoMultiplier)))
    toc
    rhoMultiplier = rhoMultiplier + 1;
end % end while loop
if numPointsRho >= 10^8
    fprintf('Rho accuracy questionable.\n')
end % end if


fprintf(successFile, 'The rotation number is approximately %.15g.\n  '...
    , rho(end));
fprintf('The rotation number is approximately %.15g.\n  ', rho(end))
rho = rho(end);
toc
fprintf('\n')



%% Initialize parameterization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param = cell(2,K);
fprintf('Compute the initial Fourier modes with %d initial modes... \n  ', initialModes)
% for k = 1:K
% 
%     
%     
%     [param{1,k}, param{2,k}] = fourierParam(tori{k}, initialModes, K, rho, perPts(:,k)); 
% end % end for loop

yorkeSanderModes = 8
for k = 1:K
   [param{1,k}, param{2,k}] = fourierParam(tori{k}, yorkeSanderModes, K, rho, perPts(:,k)); 
end % end for loop

for k = 1:K
	param{1,k} = truncate(param{1,k},initialModes);
    param{2,k} = truncate(param{2,k},initialModes);
end % for loop


toc
fprintf('\n')


% new = truncate(obj,N)
%Example:
% currentModes = initialModes;
%             for k = 1:K
%                 param{1,k} = truncate(param{1,k},currentModes);
%                 param{2,k} = truncate(param{2,k},currentModes);
%             end % for loop

%% Plot initial parameterization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(figTrajectory);
fourierCellPlot(param, 'm');
%% Flip rho if needed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rhoFlip == 1 % rhoFlip(1) == '1'
    rho = 1 - rho;
    fprintf('Flipped rho is %.15g.\n\n', rho)
else % rhoFlip == 0
    fprintf('\n')
end % end if
%% Measure the sequence space error of the initial parameterization %%%%%%%
fprintf('Compute initial parameterization error with nu = %1.1f ... \n', nu)
beta = 10^-K; % Why do we not just set it to zero? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase = evaluate(param{2,1},0);
defectError = normPhi(beta, param, alpha, rho, phase, nu);
conjError = conjugacyError(rho, param, alpha, conjMax);
fprintf('Error           = %e \n  ', defectError)
fprintf('Conjugacy Error = %e \n  ', conjError)
toc
fprintf('\n')
%% Initial Newton-like operation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Carry out newton-like method...\n')
newtonIter = 0;
while newtonIter < 10 || defectError >= 10^-errorLimit
    [betaNew, paramNew] = newtonStep(beta, param, alpha, rho, phase);
    defectErrorNew = normPhi(betaNew, paramNew, alpha, rho, phase, nu);
    if defectErrorNew < defectError% || newtonIter <= 2) && defectErrorNew < 50
        newtonIter = newtonIter + 1;
        defectError = defectErrorNew;
        beta = betaNew;
        param = paramNew;
        conjError = conjugacyError(rho, param, alpha, conjMax);
        fprintf('    Iteration %d, error = %e \n', newtonIter, defectError)
%         fprintf('    Machine Epsilon:  %e \n', eps(defectError))
        fprintf('       conjugacy error = %e \n', conjError)
    else % Reject new value
        fprintf('    Iteration %d, error limit exceeded \n',newtonIter + 1)
        if newtonIter == 0
            fprintf('  Newton failed. \n')
        end
        break
    end
end % while loop
fprintf('\n')
%% Compute initial sobolev norm and tail value %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
sobolevGrid = zeros(sobolevMax,1);
    fprintf('  Sobolev Norms:\n');
    fprintf(successFile, '  Sobolev Norms:\n');
for m = 1:sobolevMax
    maxSobolevNorm = max(sobolevNorm(param,m));
    fprintf('    k = %d norm is %e \n', m, maxSobolevNorm)
    fprintf(successFile, '    k = %d norm is %e \n', m, maxSobolevNorm);
    sobolevGrid(m, 1) = maxSobolevNorm;
end % for loop
endValue = maxEndVal(param);
fprintf('  Tail size: %e \n', endValue)
fprintf(successFile, '  Tail size: %e \n', endValue);
%% Plot initial newton output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(figTrajectory);
fourierCellPlot(param, 'k')
%% Plot initial newton output on new figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figParameterization = figure(2);
fourierCellPlot(param, 'k')
set(figParameterization, 'Units', 'Normalized', 'OuterPosition', [.25 .4 .25 .4]);
set(figParameterization, 'Name', 'Parameterizations')
%% Plot and save initial coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figCoefficients = figure(3);
set(figCoefficients, 'Units', 'Normalized', 'OuterPosition', [.47 .5 .25 .4]);
set(figCoefficients, 'Name', 'Log Plot of Coefficients');
logCoeffPlot(param)
saveas(figCoefficients,strcat('./',fname,'/initialCoeffs.png'));







