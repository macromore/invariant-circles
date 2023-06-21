clear vars
format long



N = 200000

sum = 0;
terms = [];

rho = (1+sqrt(5))/2
theta = 0.1

for n = 1:N-1

    t1 = n/N;
    t2 = (n+1)/N;
    thisWeight = exp(-1/(t1*(1-t1)));
    nextWeight = exp(-1/(t2*(1-t2)));
    thisTerm = thisWeight - nextWeight;
    
    theta_nminus = mod(theta + (n-1)*rho, 1);
    theta_n = mod(theta + n*rho, 1);
    sum = sum + (theta_n - theta_nminus)*(thisTerm);
    terms = [terms, thisTerm];
end

sum






