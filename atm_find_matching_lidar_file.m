function atm_lidar_match = atm_find_matching_lidar_file(f_name_cambot,lidar_dir,data_product_id)
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% SUMMARY:          The atm_find_matching_lidar_file function searches for the lidar waveform file that overlaps with the ATM CAMBOT image based on time tags
%                   in the files names, the data itself and bounding boxes for the CAMBOT and lidar data files.
%                   The IOCAM1B data product and documentation are available from NSIDC at https://nsidc.org/data/iocam1b.
%                   ATM laser products that can be used with this function are:
%                   https://nsidc.org/data/ILATMW1B and
%                   https://nsidc.org/data/ILNSAW1B
%               
%                   NOTE: when modifying this function for use with the Digital Mapping System (DMS) L1B data products the DMS timg tags, which are in GPS
%                         time, need to be converted to UTC time. The current version for this function does not support the use of DMS data products.
%                   
% SYNTAX:           atm_lidar_match = atm_find_matching_lidar_file(f_name_cambot,lidar_dir);
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% INPUT PARAMETERS
% 
%  f_name_cambot    File name of ATM CAMBOT L1B geotiff file to be read including full pathname.
%                   The data product that should be used with this function is IOCAM1B and is available at:
%                   https://nsidc.org/data/iocam1b
%
%  lidar_dir        Full pathname of directory with ATM lidar waveform data products.
%
%  data_product_id  Only 'ILATMW1B' and 'ILNSAW1B' are supported.
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% OUTPUT PARAMETERS
%
%  atm_lidar_match  Varible name of a MATLAB structured array (struct) containing search results.
%                   Metadata and processing information are included as well. 
%                   atm_lidar_match must be a valid MATLAB variable name.
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Author:           Michael Studinger, NASA Goddard Space Flight Center, Greenbelt MD, USA.
% Version:          1.02 - August 10, 2021
% See also:         
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------

%% parse input and output parameters and do basic error checking

% check if user has a license for the mapping toolbox 
if (license('test','map_toolbox') == 0)                      
        error('atm_find_matching_lidar_file:mapping_tool_box_license', '\n\tNo license for mapping toolbox found. Script aborted.')
end

% check number of input arguments
if nargin ~= 3   
        error('atm_find_matching_lidar_file:nargin', ['\n\tERROR:  Number of input arguments must be 2:\n', ...
        '\tSYNTAX: atm_lidar_match = atm_find_matching_lidar_file(f_name_cambot,lidar_dir);\n'])
end

% check number of output arguments
if nargout ~= 1  
        error('atm_find_matching_lidar_file:nargout', ['\n\tERROR:  Number of output arguments must be 1:\n', ...
        '\tSYNTAX: atm_lidar_match = atm_find_matching_lidar_file(f_name_cambot,lidar_dir);\n'])
end

% check if CAMBOT file exists
if (ischar(f_name_cambot) && exist(f_name_cambot) == 0)
        warndlg({'CAMBOT Input file not found!';'Script aborted.'},'!! Warning !!');
        error('atm_find_matching_lidar_file:file_chk', '\n\tCAMBOT input file not found. Script aborted.')
end

% check if lidar directory exists and is a file folder
if (ischar(lidar_dir) && exist(lidar_dir) ~= 7)
        warndlg({'Lidar directory not found!';'Script aborted.'},'!! Warning !!');
        error('atm_find_matching_lidar_file:file_chk', '\n\tLidar directory not found. Script aborted.')
end

% quietly fix missing platform-specific file separator in lidar_dir in case it's needed
if ~(strcmp(lidar_dir(end),filesep))
        lidar_dir = [lidar_dir filesep];
end

% check lidar data product identifier using case insensitve string comparison
if ~(ischar(data_product_id) && strcmpi(data_product_id,'ILATMW1B') || ischar(data_product_id) && strcmpi(data_product_id,'ILNSAW1B'))  
        error('atm_find_matching_lidar_file:data_product_id', '\n\tInvalid data product ID. Only ILATMW1B and ILNSAW1B are supported. Script aborted.')
end

%% read metadata and UTC time tag from CAMBOT natural color geotiff image
%
%  The CAMBOT L1B files are named according to the following convention:
%           1         2         3         4   
%  1234567890123456789012345678901234567890
%  IOCAM1B_YYYY_LO_NASA_yyyyddmo_HHMMSS.DDDD.tif 

[filepath,cambot_name,ext] = fileparts(f_name_cambot); clear filepath ext;
date_num_cambot = datenum(cambot_name(22:36),'yyyymmdd-HHMMSS');

% get bounding box for spatial search
info_tiff = geotiffinfo(f_name_cambot); % read parameters needed for geolocation
cambot_bbox_lat = [info_tiff.CornerCoords.Lat];
cambot_bbox_lon = wrapTo360([info_tiff.CornerCoords.Lon]);
[cambot_bbox_lat, cambot_bbox_lon] = closePolygonParts(cambot_bbox_lat, cambot_bbox_lon,'degrees');

%% search for matching ATM lidar file

%  ILATMW1B and ILNSAW1B HDF5 files are named according to the following convention:
%           1         2         3
%  12345678901234567890123456789012345
%  ILATMW1B_YYYYMMDD_hhmmss.atmNXTn.h5 
%  ILNSAW1B_YYYYMMDD_HHMMSS.atm6BT7.h5 

file_list_lidar = dir([lidar_dir data_product_id '_????????_??????.atm????.h5']);  % make sure to only get NSIDC data products

first_shot_numdate = zeros(length(file_list_lidar),1); % to speed up loop execution
last_shot_numdate  = zeros(length(file_list_lidar),1);
result_flag        = zeros(length(file_list_lidar),1); % set to 1 if CAMBOT time tag is inside window of lidar data collection

for n = 1:length(file_list_lidar)

    f_name_inp = [file_list_lidar(n).folder filesep file_list_lidar(n).name];
    
    % time of first and last shot
    first_shot_str = h5read(f_name_inp,'/ancillary_data/time/first_shot');         % 2017-05-10T13:26:28-0000 in ISO-8601 standard
    last_shot_str  = h5read(f_name_inp,'/ancillary_data/time/last_shot');
    
    first_shot_numdate(n) = datenum(first_shot_str{1},'yyyy-mm-ddTHH:MM:SS-0000'); % just in case fields are duplicated as they are sometimes.
    last_shot_numdate(n)  = datenum(last_shot_str{1}, 'yyyy-mm-ddTHH:MM:SS-0000');
    
    clear f_name_inp first_shot_str last_shot_str; % clean up memory
    
    % set flag to 1 if CAMBOT time tag is within the first and last shot of lidar data collection
    if (date_num_cambot >= first_shot_numdate(n) && date_num_cambot <= last_shot_numdate(n))
        result_flag(n) = 1;
    end
    
end

% double check to make sure this is right file using the bounding boxes of CAMBOT and lidar data. 

n_results = find(result_flag == 1);

if numel(n_results) == 1
    
    f_name_lidar = [file_list_lidar(n_results).folder filesep file_list_lidar(n_results).name];
    
    % make bounding box for lidar data
    lidar_max_lat = h5read(f_name_lidar,'/ancillary_data/meta_data/max_latitude');
    lidar_min_lat = h5read(f_name_lidar,'/ancillary_data/meta_data/min_latitude');
    lidar_max_lon = h5read(f_name_lidar,'/ancillary_data/meta_data/max_longitude');
    lidar_min_lon = h5read(f_name_lidar,'/ancillary_data/meta_data/min_longitude');
    
    lidar_bbox_lat = [lidar_min_lat lidar_max_lat lidar_max_lat lidar_min_lat lidar_min_lat];
    lidar_bbox_lon = wrapTo360([lidar_min_lon lidar_min_lon lidar_max_lon lidar_max_lon lidar_min_lon]);
    [lidar_bbox_lat, lidar_bbox_lon] = closePolygonParts(lidar_bbox_lat, lidar_bbox_lon,'degrees');
    
    % now determine if CAMBOT bounding box is inside lidar bounding box    
    [inside_lidar] = inpolygon(cambot_bbox_lon,cambot_bbox_lat,lidar_bbox_lon,lidar_bbox_lat);
    
    if all(inside_lidar == 0) % CAMBOT bounding box does not seem to overlap with lidar bounding box
        error('atm_find_matching_lidar_file:data_product_id', '\n\tNo CAMBOT coordinates are located inside the lidar bounding box. Script aborted.')
    else
        atm_lidar_match.inside_lidar = inside_lidar;
    end
    
    atm_lidar_match.f_name_lidar           = f_name_lidar;
    atm_lidar_match.utc_time_cambot        = datestr(date_num_cambot,'yyyy-mm-dd HH:MM:SS.FFF');
    atm_lidar_match.utc_time_first_shot    = datestr(first_shot_numdate(n_results),'yyyy-mm-dd HH:MM:SS.FFF');
    atm_lidar_match.utc_time_last_shot     = datestr(last_shot_numdate(n_results), 'yyyy-mm-dd HH:MM:SS.FFF');
    atm_lidar_match.lidar_bbox_lon         = lidar_bbox_lon;
    atm_lidar_match.lidar_bbox_lat         = lidar_bbox_lat;
    atm_lidar_match.cambot_bbox_lon        = cambot_bbox_lon;
    atm_lidar_match.cambot_bbox_lat        = cambot_bbox_lat;
    atm_lidar_match.date_processed         = datestr(now);
    atm_lidar_match.m_file_used            = mfilename('fullpath');
    atm_lidar_match.user                   = getenv('UserName');
    atm_lidar_match.file_list_search       = file_list_lidar;
    atm_lidar_match.CAMBOT_bbox_tiff       = info_tiff.CornerCoords;

elseif numel(n_results) > 1
    error('atm_find_matching_lidar_file:two_many_results', '\n\tFound more than one matching lidar file. Check input data. Script aborted.')
elseif numel(n_results) == 0
    error('atm_find_matching_lidar_file:no_file_found', '\n\tNo matching lidar file found. Check input data. Script aborted.')
end

end