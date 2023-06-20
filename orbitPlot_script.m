clear variables
format long


alpha = acos(0.24)
 
P = [0.4;
     0];
 
numIterates = 5000
orbit = zeros(2, numIterates+1);
orbit(:, 1) = P;
numChanges = 0;
numCrossings = 0;
sum = 0;
sum_weighted = 0;

tic

AN = 0;

for n = 0:numIterates-1
    
    t = n/numIterates;
    thisWeight = exp(1/(t*(t-1)));
    AN = AN + thisWeight;
end


for n = 1:numIterates

    oldTheta = atan2(P(2), P(1))/(2*pi);
      
    P = henon(P,alpha);
    orbit(:, n+1) = P;

    theta = atan2(P(2), P(1))/(2*pi);
    dTheta = mod(theta - oldTheta, 1);
    
    sum = sum + dTheta;
    
    t = (n-1)/numIterates;
    thisWeight = exp(1/(t*(t-1)));
    sum_weighted = sum_weighted + thisWeight*dTheta;
     
end

rho = sum/numIterates

rho_weighted = sum_weighted/AN


time = toc


pk = [0;0];
k = 30
for n = 0:numIterates-1 
    
    t = n/numIterates;
    thisWeight = exp(1/(t*(t-1)));
    pk = pk + thisWeight*orbit(:, n+1)*exp(-2*pi*1i*k*n*rho_weighted);
     
end

pk = pk/AN


% 
% figure 
% hold on 
% plot(orbit(1, 2:end)', orbit(2, 2:end)', 'k*')
% plot([0, 0.39], [0, 0], 'b', 'LineWidth', 3)
% plot(0,0, 'r*')
% plot(orbit(1, 1), orbit(2,1), 'go')
% axis([-0.5 0.5 -0.5 0.5])




