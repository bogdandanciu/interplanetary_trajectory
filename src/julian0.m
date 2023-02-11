function j0 = Julian0(year, month, day)
%   This function computes the Julian day at 0 hours UT(universal time) 
%   for any of the years between 1901 and 2099.
%
%   INPUTS:
%       year  = range 1901-2099
%       month = range 1-12
%       day   = range 1-31 
%   OUTPUT: 
%       j0    = Julian day at 0 hours UT

%% Julian day at 0 hours UT
j0 = 367*year - fix(7*(year + fix((month + 9)/12))/4) ...
     + fix(275*month/9) + day + 1721013.5;
end