function oe = oe_from_sv(R,V,mu)
%   This function computes the orbital elements (oe) from the 
%   state vector (R,V).
%   Based on Algorithm 4.2 from Orbital mechanics for engineering students,
%   2010, by H.D. Curtis 
%
%   INPUTS: 
%       R  = position vector state (km)
%       V  = velocity vector (km/s)
%       mu = gravitational parameter (km^3/s^2)
%   OUTPUT: 
%       oe = vector of orbital elements 
%
%   VARIABLES DESCRIPTION:
%       r, v - the magnitudes of R and V
%       vr   - radial velocity component (km/s)
%       H    - the angular momentum vector (km^2/s)
%       h    - the magnitude of H (km^2/s)
%       incl - inclination of the orbit (rad)
%       N    - the node line vector (km^2/s)
%       n    - the magnitude of N
%       cp   - cross product of N and R
%       RA   - right ascension of the ascending node (rad)
%       E    - eccentricity vector
%       e    - eccentricity (magnitude of E)
%       w    - argument of perigee (rad)
%       TA   - true anomaly (rad)
%       a    - semimajor axis (km)

%% Calculate distance and speed
r = norm(R);
v = norm(V);

%% Calculate the radial velocity
vr = dot(R,V)/r;

%% Calculate the specific angular momentum 
H = cross(R,V);
h = norm(H);

%% Calculate the inclination 
inclination = acos(H(3)/h);

%% Calculate the node line N
N = cross([0 0 1],H);
n = norm(N);

eps = 1.e-10;
%% Calculate the right ascension the ascending node
if n ~= 0
    RA = acos(N(1)/n);
    if N(2) < 0
        RA = 2*pi - RA;
    end
else
    RA = 0;
end

%% Calculate the eccentricity vector and magnitude 
E = 1/mu*((v^2 - mu/r)*R - r*vr*V);
e = norm(E);

%% Calculate the argument of perigee 
if n ~= 0
    if e > eps
        w = acos(dot(N,E)/n/e);
    if E(3) < 0
        w = 2*pi - w;
    end
    else
        w = 0;
    end
else
    w = 0;
end

%% Calculate the true anomaly
if e > eps
    TA = acos(dot(E,R)/e/r);
    if vr < 0
        TA = 2*pi - TA;
    end
else
    cp = cross(N,R);
    if cp(3) >= 0
        TA = acos(dot(N,R)/n/r);
    else
        TA = 2*pi - acos(dot(N,R)/n/r);
    end
end
a = h^2/mu/(1 - e^2);

%% Orbial elements vector
oe = [h e RA inclination w TA a];
end 