% Compute the same circle with successively more modes and report the
% results.
% 
% Notes: 
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
initialModes = 25;
modeStep = 25;
maxModes = 350;
% Number of points for plotting and computing initial parameterization
numPoints = 1e3;
% Do we need to replace rho with 1-rho
rhoFlip = 1;
% Error limit on rho
rhoLimit = 16;
% Error limit in Newton
errorLimit = 14;
% Error limit in tail of truncated series
tailLimit = 14;
% Sobolev space H^m max
sobolevMax = 10;
% Conjugacy check max
conjMax = 100;
%% Example 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alpha = acos(0.24);
% initialP = [0.4,0];
%% Example 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alpha = acos(0.24);
% initialP = [0.5,0];
%% Example 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = acos(-0.95);
initialP = [0,-2.65];
% numPoints = 12000;
%% Plot color %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color = 'm';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% For general use, nothing below here needs to be changed.                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ask for a file name for the diary output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = input('Name folder for output? (can be root/folder, etc.) ','s');
mkdir(fname);
diary(strcat('./',fname,'/output.txt'))
tableFile = fopen(strcat('./', fname, '/table.txt'), 'w');
%% Header information for the file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Parameters: \n')
fprintf(['P =', repmat(' ',1,22),'[ %f, %f ]\n'], initialP(1),initialP(2))
fprintf(['alpha =', repmat(' ',1,18), '%.15g', repmat(' ',1,10),...
    'probably acos(something)\n'], alpha)
fprintf('nu =                     %.15g\n', nu)
fprintf('sobolev max space k      %d\n', sobolevMax)
fprintf('initial number of modes  %d\n', initialModes)
fprintf('Mode step                %d\n', modeStep)
fprintf('points per circle        %d\n', numPoints)
fprintf('max tail size            %.15g\n', 10^-tailLimit)
fprintf('max error                %.15g\n', 10^-errorLimit)
fprintf('Conjugacy iterations     %d\n', conjMax)
% fprintf('rho difference limit     %.15g \n', 10^-rhoLimit)
fprintf('rho flip                 %d\n\n', rhoFlip)
fprintf(tableFile, '\\begin{tabular}{c||l|c|c}\n Numer of Points &  $\\rho$ Approximation & Relative Error & Machine Epsilon \\\\ \n \\hline \n');
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
axis([-2 2 -2 2]) % Do I need to do this?
%% Ask the number of tori and separate the orbits %%%%%%%%%%%%%%%%%%%%%%%%%
K = input('How many periodic tori? ');
fprintf('\nSeparate the trajecotory... \n    '    )
trajectory = pointTrajectory(initialP, alpha, numPoints*K);
tori = trajectorySeparator(trajectory, K);
toc
fprintf('\n')
%% Find periodic points of Tori %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Finding periodic points ... \n    ')
perPts = periodicPointFinder(tori, K, alpha);
toc
fprintf('\n')
plot(perPts(1,:),perPts(2,:),'.k')
%% Find rotation number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Compute the rotation number... \n')
rhoMultiplier = 1;
diffRho = 1;
numPointsRho = 0;
rho = 0;
while (diffRho > 10^-rhoLimit && numPointsRho < 10^8) || rhoMultiplier < 11
    % Start with 100 points per tori and up the multiplier each iteration
    numPointsRho = K*1000*rhoMultiplier; 
    fprintf('Number of points per tori         %d\n', 1000*rhoMultiplier)
    trajectory2 = pointTrajectory(initialP, alpha, numPointsRho);
    tori2 = trajectorySeparator(trajectory2, K);
    rho(rhoMultiplier) = weightedBirkoffRotationNumber(tori2, K, perPts);
    fprintf('  Rho approximates as %.15g. \n',rho(rhoMultiplier))
    if rhoMultiplier > 1
        diffRho = abs(rho(rhoMultiplier) - rho(rhoMultiplier-1));
        fprintf('  Delta rho:          %e \n'...
            , abs(rho(rhoMultiplier)-rho(rhoMultiplier-1)))
    end
    fprintf('  Machine epsilon is  %e \n  ', eps(rho(rhoMultiplier)))
    fprintf(tableFile, ' %d & %.15f & %g & %g \\\\ \n', 1000*rhoMultiplier, rho(rhoMultiplier), diffRho, eps(rho(rhoMultiplier)));
    toc
    rhoMultiplier = rhoMultiplier + 1;
end % end while loop
if numPointsRho >= 10^8
    fprintf('Rho accuracy questionable.\n')
end % end if
fprintf('The rotation number is approximately %.15g.\n  ', rho(end))
rho = rho(end);
toc
fprintf('\n')
fprintf(tableFile, '\\end{tabular} %% Matlab generated table');
diary off;
fclose all;