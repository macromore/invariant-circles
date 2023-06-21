function Pn = henon(P,alpha)
% Area preserving henon map
%
% Define variables:
% P           - P (input)
% alpha       - parameter of map(input)
% Pn          - Output, next point
% 
% Dependencies: 
% N/A
    xn = P(1)*cos(alpha) - P(2)*sin(alpha) + P(1)*P(1)*sin(alpha);
    yn = P(1)*sin(alpha) + P(2)*cos(alpha) - P(1)*P(1)*cos(alpha);
    Pn = [xn; yn];
end % end henon