function [oe, r, v, jd] = planet_oe_and_sv ...
(planet_id, year, month, day, hour, minute, second)
%   This function computes the orbital elements, the state vector, the velocity
%   and the julien day for a planet. 
%   Based on Algorithm 8.1 from Orbital mechanics for engineering students,
%   2010, by H.D. Curtis 
%
%   INPUTS: 
%        planet_id =  planet identifier - 1 to 9, from Mercury to Pluto
%        year      = range 1901-2099
%        month     = range 1-12
%        day       = range 1-31 
%        hour      = range 0-23
%        minute    = range 0-60
%        second    = range 0-60 
%   OUTPUTS: 
%        oe        = vector of heliocentric elements 
%        r         = heliocentric position vector (km)
%        v         = heliocentric velocity vector (km)
%        jd        = julian day number
%
%   VARIABLES DESCRITPTION 
%       mu        - gravitaional parameter of the sun (km^3/s^2)
%       deg       - conversion between degrees to radinas
%       oe        - vector comprising the heliocentric elements of the planet
%                   [h e RA incl w TA a w_hat L M E]
%                       h - angular momentum (km^2/s)
%                       e - eccentricity
%                       RA - right ascension (deg)
%                       incl - inclination (deg)
%                       w - argument of perihelion (deg)
%                       TA - true anomaly (deg)
%                       a - semimajor axis (km)
%                       w_hat - longitude of perihelion ( = RA + w) (deg)
%                       L - mean longitude ( = w_hat + M) (deg)
%                       M - mean anomaly (deg)
%                       E - eccentric anomaly (deg)
%       j0       - Julian day number of the date at 0 hr UT
%       ut       - universal time in fractions of a day
%       jd       - julian day number of the date and time
%       J2000_oe - row vector of J2000 orbital elements from Table 9.1
%       rates    - row vector of Julian centennial rates from Table 9.1
%       t0       - Julian centuries between J2000 and jd
%                  elements - orbital elements at jd

%% Constants
global mu
deg = pi/180;

%% Calculate julian day at 0 hr and universal time
j0 = Julian0(year, month, day);
ut = (hour + minute/60 + second/3600)/24;
%Calculate the julian day of the date and time
jd = j0 + ut;

%% Obtain the data for the selected planet
[J2000_oe, rates] = planetary_ephemeris(planet_id);
t0 = (jd - 2451545)/36525;
el = J2000_oe + rates*t0;
a = el(1);
e = el(2);
h = sqrt(mu*a*(1 - e^2));

%% Reduce the angular elements within the range 0 - 360 degrees:
inclination = el(3);
RA = wrapTo360(el(4));
w_hat = wrapTo360(el(5));
L = wrapTo360(el(6));
w = wrapTo360(w_hat - RA);
M = wrapTo360((L - w_hat));

%% Calculate eccentric anomaly and true anomaly
E = kepler_equation(e, M*deg);
TA = wrapTo360(2*atan(sqrt((1 + e)/(1 - e))*tan(E/2))/deg);

%% Orbital elements vector  
oe = [h e RA inclination w TA a w_hat L M E/deg];

%% Calculate the state vector
[r, v] = sv_from_oe([h e RA*deg inclination*deg w*deg TA*deg],mu);
return
end