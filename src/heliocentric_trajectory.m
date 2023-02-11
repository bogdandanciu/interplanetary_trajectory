function ...
[planet1, planet2, trajectory] = heliocentric_trajectory(departure, arrival)
%   This function determines the space vehicle's heliocentric
%   trajectory from the departure planet up to the arrival planet
%
%   INPUTS: 
%       departure  = [planet_id, year, month, day, hour, minute, second]
%       arrival    = [planet_id, year, month, day, hour, minute, second]
%   OUTPUTS: state vector and juien day for planet 1 and 2
%       planet1    = [R_p1, V_p1, jd1]
%       planet2    = [R_p2, V_p2, jd2]
%       trajectory = [V1, V2]
%
%   VARIABLES DESCRIPTION:
%       planet_id  - planet identifier - 1 to 9, from Mercury to Pluto 
%       year       - range 1901-2099
%       month      - range 1-12
%       day        - range 1-31 
%       hour       - range 0-23
%       minute     - range 0-60
%       second     - range 0-60 
%       jd1, jd2   - julian day at departure and arrival 
%       R_p1, V_p1 - state vector of planet 1 at departure 
%       R_p2, V_p2 - state vector of planet 2 at arrival 
%       R1, V1     - state vector of space vehicle at departure
%       R2, V2     - state vector of space vehicle at arrival

%% Departure
planet_id = departure(1);
year = departure(2);
month = departure(3);
day = departure(4);
hour = departure(5);
min = departure(6);
sec = departure(7);
%State vector of planet 1
[~, R_p1, V_p1, jd1] = planet_oe_and_sv ...
(planet_id, year, month, day, hour, min, sec);

%% Arrival 
planet_id = arrival(1);
year = arrival(2);
month = arrival(3);
day = arrival(4);
hour = arrival(5);
min = arrival(6);
sec = arrival(7);
%State vector of planet 2
[~, R_p2, V_p2, jd2] = planet_oe_and_sv ...
(planet_id, year, month, day, hour, min, sec);
time_of_flight = (jd2 - jd1)*24*3600;

%% The assumption of patched conics
R1 = R_p1;
R2 = R_p2;

%% Outputs
%Space vehicle's velocity at departure and arrival assuming a prograde
%trajectory
[V1, V2] = Lambert(R1, R2, time_of_flight);
planet1 = [R_p1, V_p1, jd1];
planet2 = [R_p2, V_p2, jd2];
trajectory = [V1, V2];
end 