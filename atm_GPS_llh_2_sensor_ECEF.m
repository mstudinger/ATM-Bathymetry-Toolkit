function sensor_ecef = atm_gps_llh_2_sensor_ecef(ant_lat,ant_lon,ant_ele,p, r, h,lever_arm_sensor,wgs84)
%
% SUMMARY:            The ATM_GPS_LLH_2_SENSOR_ECEF function calculates the location of an instrument's sensor using the geodetic coordinates of the 
%                     phase center of the aircraft's GPS antenna, pitch, roll, and yaw (heading) and the instrument's lever arm measured in an
%                     aircraft-fixed cartesian coordinate system. The output is in a geocentric ECEF (Earth-centered, Earth-fixed) coordinate system 
%
% SYNTAX:             sensor_ecef = atm_gps_llh_2_sensor_ecef(ant_lat,ant_lon,ant_ele,pitch, roll, heading,lever_arm_sensor,wgs84);
%
%
% HISTORY:            translated rotation matrices form IDL code from NASA's Operation IceBridge Project Science Office (Jeremy Harbeck) to MATLAB
%                   
% --------------------------------------------------------------------------------------------------------------------------------------------------------------
%
% INPUT PARAMETERS
% 
%   ant_lat           Geodetic latitude of the phase center of the aircraft's GPS antenna. Latiudes must be in decimal degrees and between ±90° latitude.
%                     This parameter is stored in /aircraft/latitude in ATM's waveform HDF5 files available from NSIDC.
%                     Data products that can be used with this function are ILATMW1B and ILNSAW1B and are available at:
%                     https://nsidc.org/data/ilatmw1b (ILATMW1B ATM's wide scanner)
%                     https://nsidc.org/data/ILNSAW1B (ILNSAW1B ATM's narrow scanner)
%
%   ant_lon           Geodetic longitude of the phase center of the aircraft's GPS antenna. Longitudes must be in decimal degrees east and
%                     between 0° and 360° longitude. This parameter is stored in /aircraft/longitude in ATM's waveform HDF5 files available from NSIDC.
%                     Data products are ILATMW1B and ILNSAW1B.
% 
%   ant_ele           Ellipsoid height (WGS-84) of the phase center of the aircraft's GPS antenna in meters.
%                     This parameter is stored in /aircraft/antenna_height in ATM's waveform HDF5 files available from NSIDC. 
%                     Data products are ILATMW1B and ILNSAW1B.
% 
%   pitch             Aircraft pitch angle in degress from IMU data. Positive pitch = aircraft nose up. 
%                     This parameter is stored in /aircraft/pitch in ATM's waveform HDF5 files available from NSDIC (data products ILATMW1B and ILNSAW1B).
%
%   roll              Aircraft roll angle in degress from IMU data. Positive roll = starboard (right) wing down.
%                     This parameter is stored in /aircraft/pitch in ATM's waveform HDF5 files available from NSDIC (data products ILATMW1B and ILNSAW1B).
%
%   heading           Aircraft heading angle in degress true north from IMU data. Positive towards east. 
%                     This parameter is stored in /aircraft/heading in ATM's waveform HDF5 files available from NSDIC (data products ILATMW1B and ILNSAW1B).
%
%   lever_arm_sensor  3x1 displacement vector of the center of the scan mirror (origin of laser measurement) measured from the phase center
%                     of the aircraft's GPS antenna towards the scan mirror in an aircraft-fixed cartesian coordinate system. 
%                     The xyz parameters (in meters) are stored in /mounting_parameters/ in ATM's waveform HDF5 files available from NSIDC 
%                     (data products ILATMW1B and ILNSAW1B):
%                     /mounting_parameters/lever_arm_x: displacement from antenna phase center to center of scan mirror toward forward (in meters)
%                     /mounting_parameters/lever_arm_y: displacement from antenna phase center to center of scan mirror toward starboard (in meters)
%                     /mounting_parameters/lever_arm_z: displacement from antenna phase center to center of scan mirror toward down (in meters)
%                     lever_arm_sensor = [lever_arm_x; lever_arm_y; lever_arm_z;];
%
%  ellipsoid          The reference spheroid for the geodetic coordinates. Length unit should be in meters.
%                     An easy way to create the "ellipsoid" argument in MATLAB is to use the command: wgs84 = wgs84Ellipsoid('meter');
%                     Alternatively, it is possible to call the ATM_GPS_LLH_2_SENSOR_ECEF function with the build in MATLAB function, which might be slower:
%                     sensor_ecef = atm_gps_llh_2_sensor_ecef(ant_lat,ant_lon,ant_ele,pitch, roll, heading,lever_arm_sensor,wgs84Ellipsoid);
%
% --------------------------------------------------------------------------------------------------------------------------------------------------------------
%
% OUTPUT PARAMETERS
%
%   sensor_ecef       3x1 column vector with cartesian xyz coordinates (in meters) of the center of the scan mirror in a geocentric Earth-centered, 
%                     Earth-fixed coordinate system (ECEF). Conversion to geodetic coordinates is done outside this function in order to vectorize 
%                     the code for faster execution.  
%
% --------------------------------------------------------------------------------------------------------------------------------------------------------------
%
% Author:             Michael Studinger, NASA Goddard Space Flight Center, Greenbelt MD, USA.
% Version:            1.0 - December 18, 2020
% --------------------------------------------------------------------------------------------------------------------------------------------------------------

%% parse input parameters and do validation checking

% check if user has a license for the mapping toolbox. consider removing this check in future versions for faster execution 
if (license('test','map_toolbox') == 0)                      
        error('atm_gps_llh_2_sensor_ecef:mapping_tool_box_license', '\n\tNo license for mapping toolbox found. Script aborted.')
end

% check number of input arguments
if (nargin ~= 8)    
        error('atm_gps_llh_2_sensor_ecef:nargin', ['\n\tERROR: Number of input arguments must be 8:\n', ...
        '\tSYNTAX: sensor_ecef = atm_gps_llh_2_sensor_ecef(ant_lat,ant_lon,ant_ele,pitch,roll,heading,lever_arm_sensor,ellipsoid);\n'])
end

% check number of output arguments
if (nargout ~= 1)
        error('atm_gps_llh_2_sensor_ecef:nargout', ['\n\tERROR: Number of out arguments must be 1:\n', ...
        '\tSYNTAX: sensor_ecef = atm_gps_llh_2_sensor_ecef(ant_lat,ant_lon,ant_ele,pitch,roll,heading,lever_arm_sensor,ellipsoid);\n'])
end

% check antenna longitude and latitude. consider removing this check in future versions for faster execution 
if any(abs(ant_lat) > 90) 
    error('atm_gps_llh_2_sensor_ecef:ant_lat', '\n\tAntenna phase center latitude (ant_lat) must be between ±90°. Script aborted.')
end
if (any(ant_lon < 0) || any(ant_lon > 360))
    error('atm_gps_llh_2_sensor_ecef:ant_lon', '\n\tAntenna phase center longitude (ant_lon) must be between 0° and 360°. Script aborted.')
end

% make sure lever_arm_sensor is a 3x1 column vector and in double precision. consider removing this check in future versions for faster execution
if (length(lever_arm_sensor) ~=3)
    error('atm_gps_llh_2_sensor_ecef:lever_arm_sensor', '\n\tlever_arm_sensor needs to be a 3x1 column vector. Script aborted.')
else
    if (iscolumn(lever_arm_sensor) == 1)     % column vector as required. make sure it's in double precision
        lever_arm_sensor = double(lever_arm_sensor);
    elseif (iscolumn(lever_arm_sensor) == 0) % means is a row vector and needs to be transformed
        lever_arm_sensor = double(lever_arm_sensor)';
    end
end

%%

[ant_ecef(1) ant_ecef(2) ant_ecef(3)] = geodetic2ecef(wgs84,ant_lat,ant_lon,ant_ele,'degrees');

%% set up rotation matrix to align sensor coordinate system with aircraft coordinate system with the x-axis in the direction of North

cp = cosd(p); sp = sind(p); cr = cosd(r); sr = sind(r); ch = cosd(h); sh = sind(h);

% the T(heading,pitch,roll) rotation aligns the sensor’s coordinate system with the aircraft coordinate system, with the x-axis in the direction of the North

T = [ ch*cp ch*sp*sr-sh*cr ch*sp*cr+sh*sr;...
      sh*cp sh*sp*sr+ch*cr sh*sp*cr-ch*sr;...
      -1.0*sp          cp*sr          cp*cr ];

% rotate the sensor coordinates into a local North, East, Down (NED) coordinate system

st = sind(ant_lat); ct = cosd(ant_lat); sl = sind(ant_lon); cl = cosd(ant_lon);

NED_R = [ -1.0*st*cl -1.0*sl -1.0*ct*cl;...
            -1.0*st*sl      cl -1.0*ct*sl;...
                    ct       0    -1.0*st; ];
                
% transform the sensor coordinates from an NED into an Earth-centered, Earth-fixed (ECEF) coordinate system

sensor_ecef = NED_R * T * lever_arm_sensor + ant_ecef'; % sensor position in ECEF coordinates

end

