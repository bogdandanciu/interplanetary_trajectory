function [planet_id, year, month, day, hour, minute, second, ...
    calendarDateStr,universalTimeStr, a_parking] = user_inputs(a)
%   This function takes the user inputs.
%    
%   INPUTS: 
%       a                = 1 for departure planet, 2 for arrival planet 
%   OUTPUTS: 
%       planet_id        = planet identifier - 1 to 9, from Mercury to Pluto 
%       year             = range 1901-2099
%       month            = range 1-12
%       day              = range 1-31 
%       hour             = range 0-23
%       minute           = range 0-60
%       second           = range 0-60 
%       calendarDateStr  = string of the calendar date
%       universalTimeStr = string of the universal time

%% Planet id
if a == 1
    fprintf('\n  Choose the planet of deprarture \n');
elseif a == 2
    fprintf('\n  Choose the planet of arrival \n');
end
fprintf('\n  1 - Mercury');
fprintf('\n  2 - Venus');
fprintf('\n  3 - Earth');
fprintf('\n  4 - Mars');
fprintf('\n  5 - Jupiter');
fprintf('\n  6 - Saturn');
fprintf('\n  7 - Uranus');
fprintf('\n  8 - Neptune');
fprintf('\n  9 - Pluto \n');    
while(1)
    planet_id = input('? ');
    if (planet_id >= 1 && planet_id <= 9)
         break;
    end
end
%% Calendar date for departure or arrival
if a == 1
    fprintf('\nInput the calendar date for departure');
elseif a == 2
    fprintf('\nInput the calendar date for arrival');
end
while(1)
    fprintf('\n(1 <= day <= 31/ 1 <= month <= 12/ year = four digits)\n');
    fprintf('(Example: 1/12/2005)\n');
    calendarDateStr = input('? ', 's');
    cdl1 = size(calendarDateStr);
    cdi1 = strfind(calendarDateStr, '/');
    % get month, day and year
    day = str2double(calendarDateStr(1:cdi1(1)-1));
    month = str2double(calendarDateStr(cdi1(1)+1:cdi1(2)-1));
    year = str2double(calendarDateStr(cdi1(2)+1:cdl1(2)));
    if (day >= 1 && day <= 31 && month >= 1 && month <= 12 ...
        && year >= 1900 && year <= 2100)
        break;
    end
end
%% Universal time for departure or arrival
if a == 1
    fprintf('\nInput the universal time for departure');
elseif a == 2
    fprintf('\nInput the universal time for arrival');
end
while(1)
    fprintf('\n(0 <= hours <= 24: 0 <= minutes <= 60: 0 <= seconds <= 60)\n');
    fprintf('(Example: 12:00:00)\n');
    universalTimeStr = input('? ', 's');
    utl1 = size(universalTimeStr);
    uti1 = strfind(universalTimeStr, ':');
    % get hours, minute and seconds
    hour = str2double(universalTimeStr(1:uti1(1)-1));
    minute = str2double(universalTimeStr(uti1(1)+1:uti1(2)-1));
    second = str2double(universalTimeStr(uti1(2)+1:utl1(2)));
    if (hour >= 0 && hour <= 24 && minute >= 0 && minute <= 60 ...
        && second >= 0 &&second <= 60) 
        break;
    end
end
%% Altitute of the departure parking orbit and of the capture orbit
if a == 1
    fprintf('\nInput the altitude of the circular departure parking orbit\n');
elseif a == 2
    fprintf('\nInput the altitude of the circular capture orbit\n');
end
while(1)
    a_parking = input('? ');
    if a_parking >= 100
        break;
    end
    fprintf('\nPlease choose a bigger altitude\n');
end
end