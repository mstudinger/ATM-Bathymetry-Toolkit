function [lat_srf,lon_srf,slant_range_srf] = atm_lake_bottom_to_surface(sns_lat,sns_lon,sns_ele,point_azimuth,point_offnadir,lake_srf_ele)
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% SUMMARY:          The MATLAB® function atm_lake_bottom_to_surface function determines the geographic coordinates and slant rage of the intersection between the
%                   laser beam and the mean lake surface elevation for laser shots that have no lake surface elevation return. 
%
% SYNTAX:           [lat_srf,lon_srf,slant_range_srf] = atm_lake_bottom_to_surface(sns_lat,sns_lon,sns_ele,point_azimuth,point_offnadir,lake_srf_ele)
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% INPUT PARAMETERS
% 
%   sns_lat         Latitudes of the center of the scan mirror. To calculate these coordinates from HDF5 waveform files use atm_GPS_llh_2_sensor_ECEF.m.
%                   Latiudes must be in decimal degrees and between ±90° latitude.
%                   Crossing the poles at ±90° has not been tested and can produce unexpected results.
%   sns_lon         Longitudes of the center of the scan mirror. To calculate these coordinates from HDF5 waveform files use atm_GPS_llh_2_sensor_ECEF.m.
%                   Longitudes must be in decimal degrees east and between 0° and 360° longitude. 
%                   Crossing the prime meridian from 0° to 360° longitude has not been tested and can produce unexpected results. 
%   sns_ele         Elevation of the laser sensor in meters above the WGS-84 reference ellipsoid. To calculate these elevations from HDF5 waveform files
%                   use atm_GPS_llh_2_sensor_ECEF.m.
%   point_azimuth   Geodetic azimuth in degress of the laser beam transmitted from the aircraft to the surface target as true compass bearing 
%                   stored in the field '/laser/point_azimuth' in an ATM HDF5 waveform file. 
%   point_offnadir  Off-nadir pointing angle in degress of the laser beam stored in the field '/laser/point_offnadir' in an ATM HDF5 waveform file.
%                   When using atm_wvfm_reader.m to import wavform data atm_wvfm.point_azimuth and atm_wvfm.point_offnadir 
%                   contain these numbers in single precision that should be converted to double precision. 
%
%                   The number of elements in the sns_lat, sns_lon, sns_ele, point_azimuth, and point_offnadir vectors must be the same.
% 
%   lake_srf_ele    Mean elevation of the lake surface in meters above the WGS-84 reference ellipsoid.
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% OUTPUT PARAMETERS
%
%   lat_srf         Variable name of MATLAB® variable containing vector with 
%   lon_srf         Variable name of MATLAB® variable containing vector with
%   slant_range_srf Variable name of MATLAB® variable containing vector with 
% 
%                   lat_srf, lon_srf, and slant_range_srf must all be a valid MATLAB® variable names.  
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Author:           Michael Studinger, NASA Goddard Space Flight Center, Greenbelt MD, USA.
% Version:          1.1 - December 2, 2021
% See also:         atm_GPS_llh_2_sensor_ECEF.m and atm_estimate_water_depth.m
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------

%% parse input parameters and do validation checking

% check if user has a license for the MATLAB® mapping toolbox.
if (license('test','map_toolbox') == 0)                      
        error('atm_lake_bottom_to_surface:mapping_tool_box_license', '\n\tNo license for MATLAB® mapping toolbox found. Script aborted.')
else
    wgs84   = wgs84Ellipsoid('m'); % MATLAB® reference ellipsoid object for the World Geodetic System 1984 (WGS-84) reference ellipsoid with axes in meters.
    sns_lon = wrapTo360(sns_lon);
end

% check number of input arguments
if (nargin ~= 6)    
        error('atm_lake_bottom_to_surface:nargin', ['\n\tERROR: Number of input arguments must be 6:\n', ...
        '\tSYNTAX: [lat_srf,lon_srf,slant_range_srf] = atm_lake_bottom_to_surface(sns_lat,sns_lon,sns_ele,point_azimuth,point_offnadir,lake_srf_ele);\n'])
end

% check number of output arguments
if (nargout ~= 3)
        error('atm_lake_bottom_to_surface:nargout', ['\n\tERROR: Number of out arguments must be 3:\n', ...
        '\tSYNTAX: [lat_srf,lon_srf,slant_range_srf] = atm_lake_bottom_to_surface(sns_lat,sns_lon,sns_ele,point_azimuth,point_offnadir,lake_srf_ele);\n'])
end

% check sensor latitute, longitude and elevation input vectors
if (~isnumeric(sns_lat) || ~isvector(sns_lat))
        error('atm_lake_bottom_to_surface:lat_sns', ['\n\tERROR: sns_lat must be a numeric vector.\n', ...
        '\tSYNTAX: [lat_srf,lon_srf,slant_range_srf] = atm_lake_bottom_to_surface(sns_lat,sns_lon,sns_ele,point_azimuth,point_offnadir,lake_srf_ele);\n'])
end

if (~isnumeric(sns_lon) || ~isvector(sns_lon))
        error('atm_lake_bottom_to_surface:sns_lon', ['\n\tERROR: sns_lon must be a numeric vector.\n', ...
        '\tSYNTAX: [lat_srf,lon_srf,slant_range_srf] = atm_lake_bottom_to_surface(sns_lat,sns_lon,sns_ele,point_azimuth,point_offnadir,lake_srf_ele);\n'])
end

if (~isnumeric(sns_ele) || ~isvector(sns_ele))
        error('atm_lake_bottom_to_surface:sns_ele', ['\n\tERROR: sns_ele must be a numeric vector.\n', ...
        '\tSYNTAX: [lat_srf,lon_srf,slant_range_srf] = atm_lake_bottom_to_surface(sns_lat,sns_lon,sns_ele,point_azimuth,point_offnadir,lake_srf_ele);\n'])
end

% check sensor longitude and latitude input values
if any(abs(sns_lat) > 90) 
    error('atm_lake_bottom_to_surface:sns_lat', '\n\tSensor latitudes (sns_lat) must be real numbers between ±90°. Script aborted.')
end

if (any(sns_lon < 0) || any(sns_lon > 360))
    error('atm_lake_bottom_to_surface:sns_lon', '\n\tSensor longitudes (sns_lon) must be real numbers between 0° and 360°. Script aborted.')
end

% check point_azimuth and point_offnadir input vectors
if (~isnumeric(point_azimuth) || ~isvector(point_azimuth))
        error('atm_lake_bottom_to_surface:point_azimuth', ['\n\tERROR: point_azimuth must be a numeric vector.\n', ...
        '\tSYNTAX: [lat_srf,lon_srf,slant_range_srf] = atm_lake_bottom_to_surface(sns_lat,sns_lon,sns_ele,point_azimuth,point_offnadir,lake_srf_ele);\n'])
end

if (~isnumeric(point_offnadir) || ~isvector(point_offnadir))
        error('atm_lake_bottom_to_surface:point_offnadir', ['\n\tERROR: point_offnadir must be a numeric vector.\n', ...
        '\tSYNTAX: [lat_srf,lon_srf,slant_range_srf] = atm_lake_bottom_to_surface(sns_lat,sns_lon,sns_ele,point_azimuth,point_offnadir,lake_srf_ele);\n'])
end

if (any(point_offnadir < 0) || any(point_offnadir > 90)) 
    error('atm_lake_bottom_to_surface:point_offnadir', '\n\tLaser beam offnadir angles (point_offnadir) must be real numbers between 0° and 90°. Script aborted.')
end

if (any(point_azimuth < 0) || any(point_azimuth > 360))
    error('atm_lake_bottom_to_surface:point_azimuth', '\n\tLaser beam azimuth angles (point_azimuth) must be real numbers between 0° and 360°. Script aborted.')
end

% check if sns_lat, sns_lon, sns_ele, point_azimuth, and point_offnadir vectors all have the same number of numerical elements
if ~isequal(length(sns_lat), length(sns_lon), length(sns_ele), length(point_azimuth), length(point_offnadir))
        error('atm_lake_bottom_to_surface:numel', '\n\tERROR: The number of elements in the sns_lat, sns_lon, sns_ele, point_azimuth, and point_offnadir vectors must be the same. Script aborted')
end

% check mean lake surface elevation input
if (numel(lake_srf_ele) ~= 1)
    error('atm_lake_bottom_to_surface:lake_srf_ele', '\n\tERROR: The mean elevation of the lake surface (lake_srf_ele) in meters above the WGS-84 reference ellipsoid must be a single real number. Script aborted')
end

%% calculate geographic coordinates of the intersection between the laser beam and the lake surface and the slant range between the sensor and the lake surface 

[lat_srf,lon_srf,slant_range_srf] = lookAtSpheroid(sns_lat,sns_lon,sns_ele - lake_srf_ele,point_azimuth,point_offnadir,wgs84);

end