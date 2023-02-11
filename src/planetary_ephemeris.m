function [J2000_oe, rates] = planetary_ephemeris(planet_id)
%
%   This function extracts the J2000 elements of a planet from a given table
%   of orbital elemets and their centennial rates.
%
%   INPUTS: 
%       planet_id = planet identifier - 1 to 9, from Mercury to Pluto
%   OUPUTS: 
%       J2000_oe  = J2000 elements of the specified planet  
%       rates     = the centennial rate of change of the J2000 elements 
%
%   VARIABLES DESCRIPTION: 
%       J2000_elements - 9x6 matrix of the J2000 elements of the nine 
%           planets Mercury to Pluto. The columms represent:
%           a = semimajor axis (AU)
%           e = eccentricity
%           i = inclination (deg)
%           RA = right ascension of the ascending node (deg)
%           w_hat = longitude of perihelion (deg)
%           L = mean longitude (deg)
%       century_rates  - 9x6 matrix where the columns represent  
%               the rates of chnge per century of the J2000 matrix
%       J2000_oe       - vector containing the J2000 orbital elements 
%               corresponding to planet_id
%       rates          - vector containg the J2000 rates corresponding 
%                        to planer_id

%% J200 elements 
J2000_elements = ...
[0.38709893  0.20563069 7.00487  48.33167  77.45645   252.25084
 0.72333199  0.00677323 3.39471  76.68069  131.53298  181.97973
 1.00000011  0.01671022 0.00005  -11.26064 102.94719  100.46435
 1.52366231  0.09341233 1.85061  49.57854  336.04084  355.45332
 5.20336301  0.04839266 1.30530  100.55615 14.75385   34.40438
 9.53707032  0.05415060 2.48446  113.71504 92.43194   49.94432
 19.19126393 0.04716771 0.76986  74.22988  170.96424  313.23218
 30.06896348 0.00858587 1.76917  131.72169 44.97135   304.88003
 39.48168677 0.24880766 17.14175 110.30347 224.06676  238.92881];

%% Rates of change per century
century_rates = ...
[0.00000066   0.00002527  -23.51 -446.30   573.57   538101628.29
 0.00000092   -0.00004938 -2.86  -996.89   -108.80  210664136.06
 -0.00000005  -0.00003804 -46.94 -18228.25 1198.28  129597740.63
 -0.00007221  0.00011902  -25.47 -1020.19  1560.78  68905103.78
 0.00060737   -0.00012880 -4.15  1217.17   839.93   10925078.35
 -0.00301530  -0.00036762 6.11   -1591.05  -1948.89 4401052.95
 0.00152025   -0.00019150 -2.09  -1681.4   1312.56  1542547.79
 -0.00125196  0.00002514  -3.64  -151.25   -844.43  786449.21
 -0.00076912  0.00006465  11.07  -37.33    -132.25  522747.90];
J2000_oe = J2000_elements(planet_id,:);
rates = century_rates(planet_id,:);
%Calculate value in km from AU
au = 149597871;
J2000_oe(1) = J2000_oe(1)*au;
rates(1) = rates(1)*au;
%Calculate value in fractions of degrees from arcseconds
rates(3:6) = rates(3:6)/3600;
end 