%% Tori hunter script
alpha = input('Alpha value? ');
phaseSpace = figure;
hold on
numTrajectories = 500;
numPoints = 250;
color = 'k';
phasespacePlot(alpha, numTrajectories, numPoints, color);
trajectory = plot(0,0,'k.');
p = [1,0];
scale = 1;
set(phaseSpace, 'Units', 'Normalized', 'OuterPosition', ...
    [.05 .5 .25 .4]);
axis([-scale, scale, -scale, scale])
fprintf('Input (0,0) to stop.\n')
while p(1) ~= 0 || p(2) ~= 0
    p(1) = input('Input x-coordinate: ');
    p(2) = input('Input y-coordinate: ');
    fprintf('\n')
    if p(1) ~= 0 || p(2) ~= 0
        delete(trajectory);
        trajectory = trajectoryPlot(p, alpha, 10000, 'g');
    end
end