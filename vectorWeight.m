function weights = vectorWeight(numPoints)
% Compute a vector of weights for a Birkhoff sum
% 
% Define variables:
% numPoints - number of weights needed (input)
% weights   - vector of weights
% 
% Dependencies: 
% N/A
    weights = linspace(0,1,numPoints);
    weights = exp(-(weights.*(1-weights)).^-2);
    weights(1) = 0;
    weights(end) = 0;
    weights = weights ./ sum(weights);
end % end vectorWeight

