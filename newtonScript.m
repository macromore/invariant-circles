% New newton, the old one has some odd pathologies where it does not
% compute the needed data first. Thus it computes things multiple times for
% no reason. This one will attempt to compute all items once, and to hard
% wire in the concept of few initial modes.3
% For now it is just a script, can/will be translated into a function in
% the future.


%P = [0; -2.65];
%alpha = acos(-0.95);

nu = 1.1;
initialmodes = 5;
modes = 100;
numpoints = 5000;

P = [0.4, 0];
%alpha = 3*pi/2;
alpha = acos(.24);

%rho = 0.757879265870744;
%  rho = 0.186334340406075;
rho = 0.206174514865714;
% 
% P = [0.5; 0];
% alpha = acos(0.24);
% nu = 1.1;
% initialmodes = 5;
% modes = 15;
% numpoints = 1000;
% % 
% P = [.5, 0];
% alpha = .8;
% % nu = 1.1;
% initialmodes = 5;
% modes = 5;
% numpoints = 1000;

% This sets additional powers of 10 for rho computation,
% numpoints*10^rhofactor
%rhofactor = 0;

format long

fname = input('Name folder for output? (can be root/folder etc.) ','s');
mkdir(fname);

diary(strcat('./',fname,'/output.txt'))

fprintf('Parameters: \n')
fprintf('P =                      [ %f, %f ]\n', P(1),P(2))
fprintf('alpha =                  %f                      probably acos(something)\n', alpha)
fprintf('nu =                     %f\n', nu)
fprintf('initial number of modes  %d\n', initialmodes)
fprintf('modes used in newton     %d\n', modes)
fprintf('number of points         %d\n', numpoints)
%fprintf('number of points for rho %d\n', numpoints*10^rhofactor)
fprintf('number of points for rho 10^8\n')
fprintf('\n')

% Compute and plot trajectory
fprintf('Computing the trajectory... ')
tic
trajectory = pointTrajectory(P, alpha, numpoints);
toc
t = toc;
fprintf('\n')

figure
hold on

plot(trajectory(1,:),trajectory(2,:),'.r');
%savefig(strcat('./',fname,'/fig1_trajectory'));

% Ask the number of tori and separate the orbits
K = input('How many periodic tori? ');
fprintf('\n')

fprintf('Separate the trajecotory... \n')
tic
tori = trajectorySeparator(trajectory, K);
t = t + toc;
fprintf('\n')

% Find periodic points of Tori
fprintf('Finding periodic points ... ')
tic
perPts = periodicPointFinder(tori, K, alpha);
toc
t = t + toc;
fprintf('\n')
plot(perPts(1,:),perPts(2,:),'.k')
%savefig(strcat('./',fname,'/fig2_perPts'));

% Find rotation number.
% fprintf('Compute the rotation number... \n')
% tic
% if rhofactor ~= 0
%     tori2 = pointTrajectory(P, alpha, numpoints*10^rhofactor);
%     tori2 = trajectorySeparator(tori2, K);
%     rho = weightedBirkoffRotationNumber(tori2, K, perPts);
% else
%     rho = weightedBirkoffRotationNumber(tori, K, perPts);
% end 
% fprintf('\n Rho is %.15g.\n\n', rho)
% toc
% t = t + toc;
% fprintf('\n')

% Rho for acos(-.95)
% rho = 0.001147323169818;
% Rho for acos(.24)
% rho = 0.190669478955379;
% Rho for .8, 10^8 points
% rho = 0.108457999567334;

% Initialize parameterization and pad with zeros
param = cell(2,K);
fprintf('Compute the initial Fourier modes with %d initial modes... ', initialmodes)
tic
for i = 1:K
   [param{1,i}, param{2,i}] = fourierParam(tori{i}, initialmodes, K, rho, perPts(:,i)); 
   param{1,i} = truncate(param{1,i},modes);
   param{2,i} = truncate(param{2,i},modes);
end
toc
t = t + toc;
fprintf('\n')

cellPeriodicPlot(param, 'm');
%savefig(strcat('./',fname,'/fig_3initial_param'));

J = input('Do we need 1-rho for the newton-like operator (0 or 1)? ');
if J(1) == 1
    rho = 1-rho;
    fprintf('\n Flipped rho is %.15g.\n\n', rho)
else
    fprintf('\n')
end

% Measure the sequence space error of the initial parameterization
fprintf('Compute initial parameterization error with nu = %1.1f ... ', nu)
tic
beta = 10^-K;
phase = evaluate(param{2,1},0,2*pi);
error = normPhi(beta, param, alpha, rho, phase, nu);
fprintf('Error = %.15g \n\n', error)
toc
t = t + toc;

% Carry out newton-like operator
fprintf('Carry out newton-like method...')
J = input('What precision to stop at, 10^-J? (suggest about 10^-14) ');
fprintf('\n')
iteration = 0;
tic

while error > 10^-J && iteration < 10
    [beta, param] = newtonStep(beta, param, alpha, rho, phase);
    iteration = iteration + 1;
    error = normPhi(beta, param, alpha, rho, phase, nu);
    fprintf('Iteration %d, error = %.15g \n', iteration, error)
    toc
    fprintf('\n')
end
maxSobolevNorm = max(sobolevNorm(param))
t = t + toc;


fprintf('Total time elapsed %f seconds. \n', t)

cellPeriodicPlot(param, 'k')
legend('Trajectory','Initial Guess','Approximation')
savefig(strcat('./',fname,'/fig4_final_param'));

diary off