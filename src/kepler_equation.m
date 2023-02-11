function E = kepler_equation(e, M)
%   This function solves Kepler's Equation using Newton's method.
%   It ouptuts the eccentric annomaly and the mean anomaly.
%
%   INPUT: 
%       e = eccentricity 
%       M = mean anomaly (radians)
%   OUTPUT: 
%       E = eccentric anomaly (radians)


%% Eccentric anomaly 
%Error tolerance:
eps = 1.e-8;
%Give a starting value for E:
if M < pi
    E = M + e/2;
else
    E = M - e/2;
end
%Iterate until E is found under the tolerance
ratio = 1;
while abs(ratio) > eps
    ratio = (E - e*sin(E) - M)/(1 - e*cos(E));
    E = E - ratio;
end
end 