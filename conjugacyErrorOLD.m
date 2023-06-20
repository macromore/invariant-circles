function conjError = conjugacyError(rho, param, alpha, numberOfPoints)
% Compute the conjugacy error of an estimated solution across the
% cohomology equations.
    K = size(param,2);
    conjError = zeros(K,numberOfPoints-1);
    for p = 1:numberOfPoints-1
        phase = p/(numberOfPoints-1);
        for k = 1:K-1
            pMap = [evaluate(param{1,k},phase); evaluate(param{2,k},phase)];
            pMap = henon(pMap, alpha);
            pConj = [evaluate(param{1,k+1},phase); evaluate(param{2,k+1},phase)];
            conjError(k,p) = norm(pMap-pConj);
        end
        pMap = [evaluate(param{1,K},phase); evaluate(param{2,K},phase)];
        pMap = henon(pMap, alpha);
        pConj = [evaluate(param{1,1},phase+K*rho); evaluate(param{2,1},phase+K*rho)];
        conjError(K,p) = norm(pMap-pConj);
    end
    conjError = max(conjError,[],'all');
end