%% INTERPLANETARY MISSION DESIGN USING THE METHOD OF PATCHED CONICS

clc;
clear;
close all;

%% CONSTANTS

%Gravitational parameter of the sun
global mu 
mu = 1.327124e11; %km^3/s^2
%Conversion factor between degrees and radians
deg = pi/180;
%Astronomical unit 
au = 149597871; %km

%Planets' name array
planets = ['Mercury'; 'Venus  '; 'Earth  '; 'Mars   '; ...
           'Jupiter'; 'Saturn '; 'Uranus '; 'Neptune'; 'Pluto  '];

while(1)
%% USER INPUTS
 
fprintf('\n  < Interplanetary mission design program >\n');
fprintf('\n  < Inputs >\n');

%...Request departure parameters 
fprintf(' \nDeparture parameters \n');
[planet_id1, year1, month1, day1, hour1, min1, sec1,calendarDateStr1,...
           universalTimeStr1, a_parking] = user_inputs(1);
%...Request arrival parameters 
fprintf(' \nDeparture parameters \n');
[planet_id2, year2, month2, day2, hour2, min2, sec2,calendarDateStr2,...
    universalTimeStr2,r_capture] = user_inputs(2);

%% CALCULATION OF MISSION PARAMETERS

%...Departure parameters
departure = [planet_id1, year1, month1, day1, hour1, min1, sec1];
%...Arrival parameters 
arrival = [planet_id2, year2, month2, day2, hour2, min2, sec2];

%Obtain orbital elements and planets' state vectors 
[oe1, r1, v1, ~] = planet_oe_and_sv(planet_id1, year1, month1,...
    day1, hour1, min1, sec1);
[oe2, r2, v2, ~] = planet_oe_and_sv(planet_id2, year2, month2,...
    day2, hour2, min2, sec2);

[oe1_prime, r1_prime, v1_prime, jd1_prime] = planet_oe_and_sv(planet_id1,...
    year2, month2, day2, hour2, min2, sec2);
[oe2_prime, r2_prime, v2_prime, jd2_prime] = planet_oe_and_sv(planet_id2,...
    year1, month1, day1, hour1, min1, sec1);

%...Interplanetary trajectory
[planet1, planet2, trajectory] = heliocentric_trajectory(departure, arrival);
%Planet1 state vector
R1 = planet1(1,1:3);
%Planet1 velocity vector
Vp1 = planet1(1,4:6);
%Planet1 julian day
jd1 = planet1(1,7);

%Planet2 state vector
R2 = planet2(1,1:3);
%Planet2 velocity vector
Vp2 = planet2(1,4:6);
%Planet2 julian day
jd2 = planet2(1,7);

%Space vehicle velocity at departure
V1 = trajectory(1,1:3);
%Space vehicle velocity at arrival
V2 = trajectory(1,4:6);
%Time of flight
tof = jd2 - jd1;

%Orbital elements of the space vehicle's trajectory based on [Rp1, V1]...
oe = oe_from_sv(R1, V1, mu);
% ...and [R2, V2]
oe3 = oe_from_sv(R2, V2, mu);
%Velocitis at infinity 
vinf1 = V1 - Vp1;
if planet_id1 < planet_id2
    vinf2 = V2 - Vp2;
else
    vinf2 = Vp2 - V2;
end

%Planet1 orbit period
T1 = 2*pi/sqrt(mu)*oe1(7)^(3/2)/3600/24;
%Planet2 orbit period
T2 = 2*pi/sqrt(mu)*oe2(7)^(3/2)/3600/24;
%Transfer orbit period
T3 = 2*pi/sqrt(mu)*oe(7)^(3/2)/3600/24;

%...Planetary departure parameters 
%Planet1 astronomical data
planet1_astronomical_data = astronomical_data(planet_id1);
%Radius of planet1
r_planet1 = planet1_astronomical_data(1);
%Gravitaional parameter of planet1
mu_planet1 = planet1_astronomical_data(3);
%Radius of the circular parking orbit
rp1 = r_planet1 + a_parking;
%Speed at the periapsis of the departure parabola
vp1 = sqrt(norm(vinf1)^2 + 2*mu_planet1/rp1);
%Speed of the circular parking orbit
vC1 = sqrt(mu_planet1/rp1);
%Delta_v required for the maneuver
delta_v_departure = vp1 - vC1;
%Eccentricity of the hyperbola at departure
e_dep = 1 + rp1*norm(vinf1)^2/mu_planet1;
%Period of the circular parking orbit
T_parking = 2*pi/sqrt(mu_planet1)*rp1^(3/2)/60;

%...Planetary arrival parameters
%Planet2 astronomical data
planet2_astronomical_data = astronomical_data(planet_id2);
%Radius of planet2
rp2 = planet2_astronomical_data(1);
%Gravitaional parameter of planet2
mu_planet2 = planet2_astronomical_data(3);
%Radius of the circular capture orbit
r_p_arrival = rp2 + r_capture;
%Speed at the periapsis of the departure parabola
vp2 = sqrt(norm(vinf2)^2+2*mu_planet2/r_p_arrival);
%Speed of the circular capture orbit
vC2 = sqrt(mu_planet2/r_p_arrival);
%Delta_v required for the maneuver
delta_v_arrival = vp2 - vC2;
%Eccentricity of the hyperbola at arrival
e_arrive = 1 + r_p_arrival*norm(vinf2)^2/mu_planet2;
%Period of the circular capture orbit
T_parking2 = 2*pi/sqrt(mu_planet2)*r_p_arrival^(3/2)/60;

%Total delta_v for the mission 
delta_v_total = delta_v_departure + delta_v_arrival;

%Error message in case of absurd choices for departure and arrival times
if delta_v_departure > 60 || delta_v_arrival > 60 || tof < 0
    clc;
    try_again = ['\n\nPlease be careful with the dates of departure'...
            ' and arrival.Depending on \nthe planets of choice the'...
            'time of flight should range from several \nmonths up to'...
            'several years and even tens of years.\n *Try again*\n'];
    fprintf(try_again);
    pause(5);
else 
    break
end
end

%% OUTPUTS

fprintf('\n  < Results >\n');

fprintf('\nDeparture planet             ');
disp(planets(planet_id1,:));
fprintf('Departure calendar date      ');
disp(calendarDateStr1);
fprintf('Departure universal time     ');
disp(universalTimeStr1);
fprintf('\nDeparture julian date        %12.6f', jd1);

fprintf('\n\nArrival planet               ');
disp(planets(planet_id2,:));
fprintf('Arrival calendar date        ');
disp(calendarDateStr2);
fprintf('Arrival universal time       ');
disp(universalTimeStr2);
fprintf('\nArrival julian date          %12.6f', jd2');

fprintf('\n\nTransfer time              %12.6f  days \n ', tof);

fprintf('\n\nHeliocentric ecliptic orbital elements of the departure planet\n');
fprintf('--------------------------------------------------------------');
fprintf ('\n        sma (AU)              eccentricity          inclination (deg)         argper (deg)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', oe1(7)/au, oe1(2), oe1(4), oe1(5));
fprintf ('\n       raan (deg)          true anomaly (deg)       longper(deg)              period (days)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', oe1(3), oe1(6), oe1(8), T1);

fprintf('\n\nHeliocentric ecliptic orbital elements of the transfer orbit\n')
fprintf('prior to reaching the sphere of influence of the arrival planet\n');
fprintf('---------------------------------------------------------------');
fprintf ('\n        sma (AU)              eccentricity          inclination (deg)         argper (deg)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', oe(7)/au, oe(2), oe(4)/deg, oe(5)/deg);
fprintf ('\n       raan (deg)          true anomaly (deg)       period (min)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', oe(3)/deg, oe(6)/deg, T3);

fprintf('\n\nHeliocentric ecliptic orbital elements of the transfer orbit\n')
fprintf('prior to reaching the sphere of influence of the arrival planet\n');
fprintf('---------------------------------------------------------------');
fprintf ('\n        sma (AU)              eccentricity          inclination (deg)         argper (deg)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', oe3(7)/au, oe3(2), oe3(4)/deg, oe3(5)/deg);
fprintf ('\n       raan (deg)          true anomaly (deg)       period (min)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', oe3(3)/deg, oe3(6)/deg, T3);

fprintf('\n\nHeliocentric ecliptic orbital elements of the arrival planet\n');
fprintf('--------------------------------------------------------------');
fprintf ('\n        sma (AU)              eccentricity          inclination (deg)         argper (deg)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', oe2(7)/au, oe2(2), oe2(4), oe2(5));
fprintf ('\n       raan (deg)          true anomaly (deg)       longper(deg)              period (days)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', oe2(3), oe2(6), oe2(8), T2);

fprintf('\n\nDeparture velocity vector and magnitude\n');
fprintf('\nx-component of departure velocity                          %12.6f  km/s',V1(1));
fprintf('\ny-component of departure velocity                          %12.6f  km/s', V1(2));
fprintf('\nz-component of departure velocity                          %12.6f  km/s', V1(3));
fprintf('\ndeparture velocity magnitude                               %12.6f  km/s',norm(V1(3)));

fprintf('\n\nArrival velocity vector and magnitude\n');
fprintf('\nx-component of arrival velocity                            %12.6f  km/s', V2(1));
fprintf('\ny-component of arrival velocity                            %12.6f  km/s', V2(2));
fprintf('\nz-component of arrival velocity                            %12.6f  km/s', V2(3));
fprintf('\ndeparture velocity magnitude                               %12.6f  km/s',norm(V2));

fprintf('\n\nHyperbolic excess velocity vector and magnitude at departure\n');
fprintf('\nx-component of hyperbolic excess velocity at departure     %12.6f  km/s',vinf1(1));
fprintf('\ny-component of hyperbolic excess velocity at departure     %12.6f  km/s', vinf1(2));
fprintf('\nz-component of hyperbolic excess velocity at departure     %12.6f  km/s', vinf1(3));
fprintf('\nhyperbolic excess velocity at departure magnitude          %12.6f  km/s',norm(vinf1));

fprintf('\n\nHyperbolic excess velocity vector and magnitude at arrival\n');
fprintf('\nx-component of hyperbolic excess velocity at arrival       %12.6f  km/s', vinf2(1));
fprintf('\ny-component of hyperbolic excess velocity at arrival       %12.6f  km/s', vinf2(2));
fprintf('\nz-component of hyperbolic excess velocity at arrival       %12.6f  km/s', vinf2(3));
fprintf('\nhyperbolic excess velocity at arrival magnitude            %12.6f  km/s',norm(vinf2));

fprintf('\n\n<Planetary departure parameters>\n');

fprintf('\nAltitude of the parking orbit                                             %12.6f  km', a_parking);
fprintf('\nPeriod of the parking orbit                                               %12.6f  min', T_parking);
fprintf('\nSpeed of the space vehicle in its circular orbit                          %12.6f  km/s', vC1);
fprintf('\nRadius to periapsis of the departure hyperbola                            %12.6f  km', rp1);
fprintf('\nEccentricity of the departure hyperbola                                   %12.6f', e_dep);
fprintf('\nSpeed of the space vehicle at the periapsis of the departure hyperbola    %12.6f  km/s', vp1);
fprintf('\nDelta_v for departure                                                     %12.6f  km/s', delta_v_departure);

fprintf('\n\n<Planetary rendezvous parameters>\n');

fprintf('\nAltitude of the capture orbit                                             %12.6f  km', r_capture);
fprintf('\nPeriod of the capture orbit                                               %12.6f  min', T_parking2);
fprintf('\nSpeed of the space vehicle in its circular orbit                          %12.6f  km/s', vC2);
fprintf('\nRadius to periapsis of the arrival hyperbola                              %12.6f  km', r_p_arrival);
fprintf('\nEccentricity of the arrival   hyperbola                                   %12.6f', e_arrive);
fprintf('\nSpeed of the space vehicle at the periapsis of the arrival hyperbola      %12.6f  km/s', vp2);
fprintf('\nDelta_v for arrival                                                       %12.6f  km/s', delta_v_arrival);

fprintf('\n\nTotal delta_v for the mission                                           %12.6f  km/s\n', delta_v_total);


%% 3D GRAPHICAL REPRESENTATION OF THE HELIOCENTRIC TRAJECTORY

oea = zeros(360,6);oeb = zeros(360,6);oec = zeros(360,6);
ra = zeros(360,3);va = zeros(360,3);
rb = zeros(360,3);vb = zeros(360,3);
rc = zeros(360,3);vc = zeros(360,3);
rxa = zeros(360,1);rya = zeros(360,1);rza = zeros(360,1);
rxb = zeros(360,1);ryb = zeros(360,1);rzb = zeros(360,1);
rxc = zeros(360,1);ryc = zeros(360,1);rzc = zeros(360,1);
n=1;

 for i = 1:360
    oea(i,:)=[oe1(1), oe1(2), oe1(3)*deg, oe1(4)*deg, oe1(5)*deg,...
        oe1(6)*deg+i*deg];
    oeb(i,:)=[oe2(1), oe2(2), oe2(3)*deg, oe2(4)*deg, oe2(5)*deg,...
        oe2(6)*deg+i*deg];
    oec(i,:)=[oe(1), oe(2), oe(3), oe(4), oe(5), oe(6)+i*deg];
    
    [ra(i,:), va(i,:)] = sv_from_oe(oea(i,:), mu);
    [rb(i,:), vb(i,:)] = sv_from_oe(oeb(i,:), mu);
    [rc(i,:), vc(i,:)] = sv_from_oe(oec(i,:), mu);
    
    rxa(n) = ra(i,1);
    rya(n) = ra(i,2);
    rza(n) = ra(i,3);
    
    rxb(n) = rb(i,1);
    ryb(n) = rb(i,2);
    rzb(n) = rb(i,3);
    
    rxc(n) = rc(i,1);
    ryc(n) = rc(i,2);
    rzc(n) = rc(i,3);
    
    n = n+1;
 end

%Set color to black 
colordef black;

hold on
grid on
axis equal
%Initial orbit
plot3(rxa/au,rya/au,rza/au,'-r','LineWidth', 1.5);
%Final orbit
plot3(rxb/au,ryb/au,rzb/au,'-g','LineWidth', 1.5);
%Transfer orbit
if oe(6) > oe3(6)
    a = (floor(oe3(6)/deg)+(360-floor(oe(6)/deg)));
else
    a = floor(oe3(6)/deg - oe(6)/deg);
end
plot3(rxc(1:a)/au,ryc(1:a)/au,rzc(1:a)/au,'-b','LineWidth', 1.5);
plot3(rxc(a:end)/au,ryc(a:end)/au,rzc(a:end)/au,'--b','LineWidth', 1.5);
%Intial and final positions of the planets
plot3(r1(1)/au,r1(2)/au,r1(3)/au,'ok','MarkerSize',7,'MarkerFaceColor','c');
plot3(r1_prime(1)/au,r1_prime(2)/au,r1_prime(3)/au,'ok','MarkerSize',7,...
    'MarkerFaceColor','r');
plot3(r2_prime(1)/au,r2_prime(2)/au,r2_prime(3)/au,'ok','MarkerSize',7,...
    'MarkerFaceColor','c');
plot3(r2(1)/au,r2(2)/au,r2(3)/au,'ok','MarkerSize',7,'MarkerFaceColor','r');
%Draw sun
plot3(0,0,0,'oy', 'MarkerSize',18,'MarkerFaceColor','y');
%Plot and label vernal equinox line
    x_vernal = oe1(7)/au;
if oe1(7) < oe2(7)
    x_vernal = oe2(7)/au;
end    
line ([0.06, x_vernal], [0, 0], 'Color', 'w');
text(1.05 * x_vernal, 0, '\Upsilon');
%Label plot and axes
title('Heliocentric Trajectory', 'FontSize', 16);
xlabel('X coordinate (AU)', 'FontSize', 14);
ylabel('Y coordinate (AU)', 'FontSize', 14);
%Legend 
lgd = legend(strcat(planets(planet_id1,:),' orbit'),...
             strcat(planets(planet_id2,:),' orbit'),...
             'Transfer trajectory','Transfer orbit',...
             ['Initial position of ' planets(planet_id1,:)],...
             ['Final position of ' planets(planet_id1,:)],...
             ['Initial position of ' planets(planet_id2,:)], ...
             ['Final position of ' planets(planet_id2,:)],...
             'Sun','Vernal equinox direction');
lgd.FontSize = 12;
pos = get(lgd,'position');
set(lgd,'position',[0.65 0.5 pos(3:4)]);
%Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%Enable 3d rotation
rotate3d on;