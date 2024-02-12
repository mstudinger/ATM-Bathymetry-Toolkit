function [utc_to_gps] = atm_GPS_to_UTC(t_gps)

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% SUMMARY:          This MATLAB® function determines the varying time offset in seconds between GPS time, a continuous time scale without leap seconds,
%                   and the Coordinated Universal Time (UTC) that includes leap seconds. UTC = GPS - offset.
%                   
%                   When a new leap second is added to the UTC time the offset between GPS and UTC times changes and is published here:
%                   https://confluence.qps.nl/qinsy/en/utc-to-gps-time-correction-32245263.html#UTCtoGPSTimeCorrection-UTCtoGPSTimeConversionTable
%                   and here:
%                   https://maia.usno.navy.mil/ser7/tai-utc.dat for TAI (International Atomic Time)
%
%                   For announcements of future leap seconds see: 
%                   https://hpiers.obspm.fr/iers/bul/bulc/bulletinc.dat
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
%                   
% SYNTAX:           utc_to_gps = atm_GPS_to_UTC(t_gps);
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% INPUT PARAMETER
% 
%  t_gps            GPS-based time tag in MATLAB®'s serial date number format. MATLAB® serial date numbers represent the whole and fractional number of days
%                   since January 0, 0000 in the proleptic ISO calendar. 
%
%                   For conversions use the MATLAB® function datenum:
%                   t_gps = datenum(Y,M,D,H,MN,S);
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% OUTPUT PARAMETER
%
%  utc_to_gps       Offset in seconds between GPS and UTC time for the input GPS time tag t_gps.  
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Author:           Michael Studinger, NASA Goddard Space Flight Center, Greenbelt MD, USA.
% Version:          1.03 - March 3, 2022
% Last updated:     03/03/2022
% See also:         atm_wvfm_info.m
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------

% check number of input arguments
if nargin ~= 1   
        error('atm_GPS_to_UTC:nargin', ['\n\tERROR:  Number of input arguments must be 1:\n', ...
        '\tSYNTAX: utc_to_gps = atm_GPS_to_UTC(t_gps);\n'])
end

% check number of output arguments
if nargout ~= 1  
        error('atm_GPS_to_UTC:nargout', ['\n\tERROR:  Number of output arguments must be 1:\n', ...
        '\tSYNTAX: utc_to_gps = atm_GPS_to_UTC(t_gps);\n'])
end

% check if t_gps is a MATLAB® serial date after the introduction of leap seconds and a date and time that is not in the future
if (~isfinite(t_gps) || t_gps < datenum(1990,01,01) || t_gps > datenum(now) + 1)
        error('atm_GPS_to_UTC:t_gps', ['\n\tERROR: t_gps must be a MATLAB® serial date number after January 1, 1990 and not in the future.\n',...
            '\tSYNTAX: utc_to_gps = atm_GPS_to_UTC(t_gps);\n']);
end

%% determine GPS to UTC offset

if (t_gps >= datenum(1990,01,01) & t_gps < datenum(1991,01,01))
    utc_to_gps = 6;
elseif (t_gps >= datenum(1991,01,01) & t_gps < datenum(1992,07,01))
    utc_to_gps = 7;    
elseif (t_gps >= datenum(1992,07,01) & t_gps < datenum(1993,07,01))
    utc_to_gps = 8;
elseif (t_gps >= datenum(1993,07,01) & t_gps < datenum(1994,07,01))
    utc_to_gps = 9;      
elseif (t_gps >= datenum(1994,07,01) & t_gps < datenum(1996,01,01))
    utc_to_gps = 10;      
elseif (t_gps >= datenum(1996,01,01) & t_gps < datenum(1997,07,01))
    utc_to_gps = 11;
elseif (t_gps >= datenum(1997,07,01) & t_gps < datenum(1999,01,01))
    utc_to_gps = 12;         
elseif (t_gps >= datenum(1999,01,01) & t_gps < datenum(2006,01,01))
    utc_to_gps = 13;
elseif (t_gps >= datenum(2006,01,01) & t_gps < datenum(2009,01,01))
    utc_to_gps = 14;
elseif (t_gps >= datenum(2009,01,01) & t_gps < datenum(2012,07,01))
    utc_to_gps = 15;        
elseif (t_gps >= datenum(2012,07,01) & t_gps < datenum(2015,07,01))
    utc_to_gps = 16;        
elseif (t_gps >= datenum(2015,07,01) & t_gps < datenum(2016,12,31)) 
    utc_to_gps = 17;
elseif (t_gps >= datenum(2017,01,01) & t_gps < datenum(2022,12,31)) % update here when future leap seconds are introduced
    utc_to_gps = 18;             
else
    warndlg({'Unable to determine UTC to GPS offset.';'Please verify GPS input time or update GPS to UTC offsets in the atm_GPS_to_UTC MATLAB® function.'},'!! Error !!');
    error('atm_GPS_to_UTC:gps_to_utc_offset', '\n\tUnable to determine UTC to GPS offset. Please verify GPS input time or update GPS to UTC offsets in teh atm_GPS_to_UTC MATLAB® function. Script aborted.')
end

end


    