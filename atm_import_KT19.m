function [dgps_utc_serial_date,srf_temp,int_temp,emissivity,response_time] = atm_import_KT19(f_name_kt19)
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% SUMMARY:          The atm_import_kt19 function imports ATM KT-19 ASCII data files available from NSIDC.
%                   The ATM KT-19 data product contains surface temperature measurements acquired by the Heitronics KT19.85 Series II Infrared Radiation Pyrometer.
%                   The ATM KT-19 data product and documentation are available from NSIDC at https://nsidc.org/data/IAKST1B/versions/2.
%                   The ATM KT-19 data product that can be used with this function is:
%                   https://nsidc.org/data/IAKST1B
%                   
% SYNTAX:           [dgps_utc_serial_date,srf_temp,int_temp,emissivity,response_time] = atm_import_kt19(f_name_kt19);
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% INPUT PARAMETER
% 
%  f_name_kt19      File name of ATM KT-19 ASCII data file to be read including full pathname.
%                   The ATM data product that should be used with this function is IAKST1B which is available at:
%                   https://nsidc.org/data/IAKST1B
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% OUTPUT PARAMETER
%
%  dgps_utc_serial_date: array of time tags for KT-19 measurements in Coordinated Universal Time (UTC) as MATLAB® serial date numbers. MATLAB® serial date numbers
%                   represent the whole and fractional number of days since January 0, 0000 in the proleptic ISO calendar.
%
%                   For all other output paramters see IAKST1B user guide at NSIDC at:
%                   https://nsidc.org/data/IAKST1B
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Author:           Michael Studinger, NASA Goddard Space Flight Center, Greenbelt MD, USA.
% Version:          1.02 - August 10, 2021
% See also:         
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% KT-19 ASCII file format:
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% # KT19 Processed Output File
% # 
% # Acquisition Start (UTC, machine time): 2019/05/12 11:46:38
% # 
% # KT19 Settings:
% # 	Emissivity constant: 0.97
% # 	Temperature Units: C
% # 	Response Time: 0.3 sec
% # 	Ambient Measurement: Internal Reference
% # 
% # Year,Day_Of_Year, Seconds_Of_Day(UTC), Latitude(Deg), Longitude(Deg), Aircraft_Altitude_Above_Ellipsoid(m), KT19_Temperature(C), KT19_Internal_Temperature(C)
% 2019,132,42400.65,67.012596,309.302641,82.897,13.56,14.16
% 2019,132,42400.75,67.012593,309.302640,82.904,13.79,14.16
% 2019,132,42400.85,67.012591,309.302640,82.911,13.85,14.16
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------

% check number of input arguments
if nargin ~= 1   
        error('atm_import_KT19:nargin', ['\n\tERROR:  Number of input arguments must be 1:\n', ...
        '\tSYNTAX: [dgps_utc_serial_date,srf_temp,int_temp,emissivity,response_time] = atm_import_KT19(f_name_kt19);\n'])
end

% check number of output arguments
if nargout ~= 5  
        error('atm_import_KT19:nargout', ['\n\tERROR:  Number of output arguments must be 5:\n', ...
        '\tSYNTAX: [dgps_utc_serial_date,srf_temp,int_temp,emissivity,response_time] = atm_import_KT19(f_name_kt19);\n'])
end

% check if KT-19 input file exists
if (ischar(f_name_kt19) && exist(f_name_kt19) == 0)
        warndlg({'KT-19 input file not found!';'Script aborted.'},'!! Error !!');
        error('atm_import_KT19:file_chk', '\n\tKT-19 input file not found. Script aborted.')
end

%%

fid = fopen(f_name_kt19);

counter = 0;

while fid
    
    counter = counter + 1;
    
    t_line = fgetl(fid);
    
    if (t_line(1) == '#') % need to modify comment lines starting with a semicolon
        
        if (counter == 6)
            emissivity = cell2mat(textscan(char(t_line(25:end)),'%f','Delimiter',' '));
        elseif (counter == 8)
            response_time = cell2mat(textscan(char(t_line(19:end)),'%f','Delimiter',' '));    
        end

    else
        break
    end
end

n_headerlines = 11; % including variable names

opts = detectImportOptions(f_name_kt19);
opts.PreserveVariableNames = 1;

T = readtable(f_name_kt19,opts);

year = table2array(T(:,1));
day  = table2array(T(:,2));
sec  = table2array(T(:,3));
lat  = table2array(T(:,4));
lon  = table2array(T(:,5));
srf_temp = table2array(T(:,7));
int_temp = table2array(T(:,8));

lon = wrapTo180(lon);

date_num_yr = datenum(year, 1, 1);
date_num_day = date_num_yr - 1 + day;
dgps_utc_serial_date = date_num_day + sec/86400;

% this code can be used to convert MATLAB serial date numbers to human readable character strings in the 'yyyy/mm/dd HH:MM:SS.FFF' format if desired
% dgps_time_stamp = datestr(dgps_utc_serial_date,'yyyy/mm/dd HH:MM:SS.FFF');

end

