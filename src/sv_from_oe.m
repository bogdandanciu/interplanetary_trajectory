function [r, v] = sv_from_oe(oe,mu)
%   This function computes the state vector using the
%   classical orbital elements (oe).
%   Based on Algorithm 4.5 from Orbital mechanics for engineering students,
%   2010, by H.D. Curtis 
%
%   INPUTS:
%       oe = vector of orbital elements 
%       mu = gravitational parameter (km^3/s^2)
%   OUTPUTS:
%       r  = position vector (km)
%       v  = velocity vector (km/s)
%
%   VARIABLES DESCRIPTION:
%       oe   - orbital elements [h e RA incl w TA] where:
%           h    - angular momentum (km^2/s)
%           e    - eccentricity
%           RA   - right ascension of the ascending node (rad)
%           incl - inclination of the orbit (rad)
%           w    - argument of perigee (rad)
%           TA   - true anomaly (rad)
%       R3_w - Rotation matrix about the z-axis through the angle w
%       R1_i - Rotation matrix about the x-axis through the angle i
%       R3_W - Rotation matrix about the z-axis through the angle RA
%       Q_pX - Matrix of the transformation from perifocal to geocentric
%              equatorial frame
%       rp   - position vector in the perifocal frame (km)
%       vp   - velocity vector in the perifocal frame (km/s)
%       r    - position vector in the geocentric equatorial frame (km)
%       v    - velocity vector in the geocentric equatorial frame (km/s)       

%% Obtain orbital elements 
h = oe(1);
e = oe(2);
RA = oe(3);
incl = oe(4);
w = oe(5);
TA = oe(6);

%% Calculate position and velocity vector in the perifocal frame
rp = (h^2/mu) * (1/(1 + e*cos(TA))) * (cos(TA)*[1;0;0] + sin(TA)*[0;1;0]);
vp = (mu/h) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);

%% Rotation matrices
R3_W = [ cos(RA) sin(RA) 0
        -sin(RA) cos(RA) 0
        0        0       1];
R1_i = [1 0         0
        0 cos(incl) sin(incl)
        0 -sin(incl) cos(incl)];
R3_w = [ cos(w) sin(w) 0
        -sin(w) cos(w) 0
        0       0      1];

%% Transformation matrix
Q_pX = (R3_w*R1_i*R3_W)';

%% Caculate position and velocity vectors in geocentric equatorial frame
r = Q_pX*rp;
v = Q_pX*vp;
%Convert the state vector (r,v) into row vectors:
r = r';
v = v';
end