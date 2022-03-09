function [atm_depth] = atm_estimate_water_depth(lon_s,lat_s,ele_s,beam_off_nadir,beam_azimuth,delta_t,verbose)
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% SUMMARY:          The ATM_ESTIMATE_WATER_DEPTH function calculates geographic locations and elevations for the returns from water-ice interfaces.
%                   The locations and elevations are calculated using an implementation of Snell's law in geographic coordinates that accounts for the 
%                   refractive index of water in both, geolocation and slant range.
%
% SYNTAX:           [atm_depth] = atm_estimate_water_depth(lon_s,lat_s,ele_s,beam_off_nadir,beam_azimuth,delta_t,verbose);
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% INPUT PARAMETERS
% 
%   lon_s           Array with longitudes of the water surface elevations.
%                   Longitudes must be in decimal degrees and between 0° and 360° or ±180° degrees. 
%                   Crossing the prime meridian from 0° to 360° longitude has not been tested and can produce unexpected results.
%
%   lat_s           Array with longitudes of the water surface elevations.
%                   Latiudes must be in decimal degrees and between ±90° latitude.
%                   Crossing the poles at ±90° has not been tested and can produce unexpected results.
%
%   ele_s           Array with elevations of the water surface.
%                   Elevations must be in meters above WGS-84 reference ellipsoid.
%
%   beam_off_nadir  Off-nadir pointing angle of laser beam. 
%                   Off-nadir angles must be in decimal degrees and between 0° and 90° degrees.
%                   This parameter is stored in /laser/point_offnadir in ATM's waveform HDF5 files available from NSIDC.
%                   Data products that can be used with this function are ILATMW1B and ILNSAW1B and are available at:
%                   https://nsidc.org/data/ilatmw1b (ILATMW1B ATM's wide scanner)
%                   https://nsidc.org/data/ILNSAW1B (ILNSAW1B ATM's narrow scanner)
%
%   beam_azimuth    Geodetic azimuth of laser beam transmitted from the aircraft to the surface target as true compass bearing.
%                   Azimuth angles must be in decimal degrees and between 0° and 360° or ±180° degrees. 
%                   Crossing the prime meridian from 0° to 360° longitude has not been tested and can produce unexpected results.
%                   This parameter is stored in /laser/point_azimuth in ATM's waveform HDF5 files available from NSIDC. 
%                   Data products are ILATMW1B and ILNSAW1B.
%
%   delta_t         Two-way travel time difference between the returns from the air-water and water-ice interfaces.
%                   Two-way travel time differences must be in nano seconds. 
%
%                   All input arrays must have the same number of elements.
% 
%   verbose         Must be logic 0 or 1. 
%                   1 displays some basic results on the console. 0 is silent.
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% OUTPUT PARAMETERS
%
%   atm_depth       Varible name of MATLAB structured array (struct) containing results.
%                   Metadata and processing information are included as well.
%                   atm_depth must be a valid MATLAB variable name.  
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Author:           Michael Studinger, NASA Goddard Space Flight Center, Greenbelt MD, USA.
% Version:          1.02 - August 6, 2021
% See also:         
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------

tStart = tic;       % start stopwatch timer

% define some processing parameters
n_water = 1.333;          % refractive index commonly used for fresh water.
n_air   = 1.0003;         % refractive index commonly used for air at STP (https://en.wikipedia.org/wiki/List_of_refractive_indices)
n_vac   = 1.0;            % refractive index in vaccum.
v_vac   = 299792458*1E-9; % speed of light in vacuum in meters per nano second

% use Snell's law to calculate speed of light in air
v_air = (n_vac/n_air)*v_vac; % speed of light in air in meters per nano seconds

%% check input data and MATLAB license

% check if user has a license for the mapping toolbox 
if (license('test','map_toolbox') == 0)                      
        error('ATM_ESTIMATE_WATER_DEPTH:mapping_tool_box_license', '\n\tNo license for mapping toolbox found. Script aborted.')
end

% check number of input arguments
if (nargin ~= 7)    
        error('ATM_ESTIMATE_WATER_DEPTH:nargin', ['\n\tERROR: Number of input arguments must be 7:\n', ...
        '\tSYNTAX: [atm_depth] = atm_estimate_water_depth(lat_s,lon_s,ele_s,beam_off_nadir,beam_azimuth,delta_t,verbose);\n'])
end

% check number of output arguments
if (nargout ~= 1)
        error('ATM_ESTIMATE_WATER_DEPTH:nargout', ['\n\tERROR: Number of out arguments must be 1:\n', ...
        '\tSYNTAX: [atm_depth] = atm_estimate_water_depth(lat_s,lon_s,ele_s,beam_off_nadir,beam_azimuth,delta_t,verbose);;\n'])
end

% determine if console output is desired
if (length(nargin) == 7)
    verbose = nargin(7);
    if (verbose ~= 0 && verbose ~= 1)
        error('ATM_ESTIMATE_WATER_DEPTH:verbose', '\n\tERROR: verbose must be either 0 or 1.')
    end
elseif(nargin == 1 || exist('verbose') == 0) % verbose has not been set 
    verbose = 0;
end

% check length of input arrays to make sure they are all the same
if (any([length(lon_s) length(lat_s) length(ele_s) length(beam_off_nadir) length(beam_azimuth) length(delta_t)] ~= length(lon_s)))
    error('ATM_ESTIMATE_WATER_DEPTH:input_arrays', '\n\tERROR: input arrays must all have the same length.')
end

% check latitude array
if (any(~isfinite(lat_s)) || ~isvector(lat_s) || any(abs(lat_s) > 90))
    error('ATM_ESTIMATE_WATER_DEPTH:lat_s', '\n\tERROR: lat_s values must be finite and between ±90°.\n');
end

% check longitude array
lon_s = wrapTo360(lon_s);
if (any(lon_s < 0) || any(lon_s > 360) || any(~isfinite(lon_s)) || ~isvector(lon_s))
    error('ATM_ESTIMATE_WATER_DEPTH:lon_s', '\n\tERROR: lon_s values must be finite and between 0° and 360° or ±180°.\n');
end

% check ele_s array
if (any(~isfinite(ele_s)) || ~isvector(ele_s))
    error('ATM_ESTIMATE_WATER_DEPTH:ele_s', '\n\tERROR: ele_s values must be finite and in meters.\n');
end

% check beam_off_nadir array - since this is coming from the HDF5 file it is unlikely to have faulty elements
if (any(~isfinite(beam_off_nadir)) || ~isvector(beam_off_nadir) || any(beam_off_nadir < 0) || any(beam_off_nadir >90))
    error('ATM_ESTIMATE_WATER_DEPTH:beam_off_nadir', '\n\tERROR: beam_off_nadir values must be finite and between 0° and 90°.\n');
end

% check beam_azimuth array - since this is coming from the HDF5 file it is unlikely to have faulty elements
beam_azimuth = wrapTo360(beam_azimuth);
if (any(~isfinite(beam_azimuth)) || ~isvector(beam_azimuth) || any(beam_azimuth < 0) || any(beam_azimuth > 360))
    error('ATM_ESTIMATE_WATER_DEPTH:beam_azimuth', '\n\tERROR: beam_azimuth values must be finite and between 0° and 360° or ±180°.\n');
end

% check delta_t array - most likely to contain NaNs
if (any(~isfinite(delta_t)) || ~isvector(delta_t))
    error('ATM_ESTIMATE_WATER_DEPTH:delta_s', '\n\tERROR: delta_t values must be finite and in nano seconds.\n');
end

%% calculate lake depth

% use Snell's law to calculate speed of light in water
v_water = (n_air/n_water)*v_air; % speed of light in water in meters per nano seconds

% use Snell's law to calculate angle of refraction in water in water
refracted_angle_w = rad2deg(asin(sin(deg2rad(beam_off_nadir) * n_air/n_water))); % refracted angle in water measured from normal

% calculate slant range in water from two-way travel time difference
slant_range_w = 0.5 * delta_t * v_water; % in meters

% calculate locations and elevations of bottom returns
% elevation angle needed as input for MATLAB function aer2geodetic is the angle in degrees from horizontal, positive up
ele_angle = -double(90 - refracted_angle_w);

% azimuth-elevation-range (AER) 
[lat_b,lon_b,ele_b] = aer2geodetic(beam_azimuth,ele_angle,slant_range_w,lat_s,lon_s,ele_s,wgs84Ellipsoid);

% % manual calculation for verification
% %water_depth = slant_range * cos(deg2rad(theta2));
% % water_depth = slant_range * cos(deg2rad(test_angle_water));

water_depth = ele_s - ele_b; % positive

% populate the waveform struct - this is for some reason not much faster than doing it inside the loop
atm_depth = struct('lon_b',lon_b,...
    'lat_b',lat_b,... 
    'ele_b',ele_b,...
    'refracted_angle_w',refracted_angle_w,...
    'slant_range_w',slant_range_w,...
    'water_depth',water_depth);   

tStart = tic;

t_elapsed = toc(tStart);
 
%% create info field containing basic metadata about the processing
atm_depth.info.date_processed         = datestr(now);
atm_depth.info.m_file_used            = mfilename('fullpath');
atm_depth.info.user                   = getenv('UserName');

%% console output if selected
if (verbose == 1)
    
    fprintf('\n-----------------------------------------------------------------------\n');
    fprintf('Processing results\n');
    fprintf('-----------------------------------------------------------------------\n');
    fprintf('Processing time (sec):               %5.2f\n',t_elapsed);
    fprintf('Minimum water depth (m):              %.2f\n',min(water_depth));
    fprintf('Maximum water depth (m):              %.2f\n',max(water_depth));
    fprintf('Mean    water depth (m):              %.2f\n',mean(water_depth));
    fprintf('Number of laser shots processed:  %8d\n',length(lon_s));
    fprintf('-----------------------------------------------------------------------\n');
    
end

end