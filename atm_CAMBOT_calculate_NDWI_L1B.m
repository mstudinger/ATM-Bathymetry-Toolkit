function varargout = atm_CAMBOT_calc_NDWI_L1B(f_name_cambot,varargin)
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% SUMMARY:          The atm_CAMBOT_calc_NDWI_L1B function calculates the normalized difference water index (NDWI), which has been modified for ice (NDWI_ice)
%                   using ATM's geolocated and orthorectified natural color (RGB) L1B data product as input. The IOCAM1B data product and documentation 
%                   are available from NSIDC at https://nsidc.org/data/iocam1b.
%
%                   The calculation of NDWI_ice is described in 
%                   Yang, K., and Smith, L. C.: Supraglacial Streams on the Greenland Ice Sheet Delineated From Combined Spectral–Shape Information 
%                   in High-Resolution Satellite Imagery, IEEE Geosci. Remote Sens. Lett., 10, 801-805, 10.1109/LGRS.2012.2224316, 2013.
%
% SYNTAX:           [ndwi,f_name_out] = atm_CAMBOT_calc_NDWI_L1B(f_name_cambot,save_geotiff); if geotiff output of NDWI_ice is desired
%                   [ndwi]            = atm_CAMBOT_calc_NDWI_L1B(f_name_cambot);              no geotiff output
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% INPUT PARAMETERS
% 
%   f_name_cambot   File name of ATM CAMBOT L1B geotiff file to be read including full pathname.
%                   The data product that should be used with this function is IOCAM1B and is available at:
%                   https://nsidc.org/data/iocam1b
%
%   save_geotiff    Flag for exporting NDWI_ice geotiff file. Must be logic 0 or 1. 
%                   1 exports an NDWI_ice geotiff file. 0 does not. If only f_name_cambot is provided as input argument no geotiff output is saved. 
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% OUTPUT PARAMETERS
%
%   ndwi            NDWI_ice matrix with the same dimensions as the geotiff natural color input image.
%                   ndwi must be a valid MATLAB variable name.
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Author:           Michael Studinger, NASA Goddard Space Flight Center, Greenbelt MD, USA.
% Version:          1.02 - July 29, 2021
% See also:         
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------

%% parse input and output parameters and do basic error checking

% check if user has a license for the mapping toolbox 
if (license('test','map_toolbox') == 0)                      
        error('atm_CAMBOT_calc_NDWI_L1B:mapping_tool_box_license', '\n\tNo license for mapping toolbox found. Script aborted.')
end

% check number of input arguments
if ~any(nargin == [1 2])   
        error('atm_CAMBOT_calc_NDWI_L1B:nargin', ['\n\tERROR:  Number of input arguments must be 1 or 2:\n', ...
        '\tSYNTAX: [ndwi,f_name_out] = atm_CAMBOT_calc_NDWI_L1B(f_name_cambot,save_geotiff);\n',...
        '\tSYNTAX: [ndwi]            = atm_CAMBOT_calc_NDWI_L1B(f_name_cambot);\n'])
end

% check number of output arguments
if ~any(nargout == [1 2])  
        error('atm_CAMBOT_calc_NDWI_L1B:nargout', ['\n\tERROR:  Number of output arguments must be 1 or 2:\n', ...
        '\tSYNTAX: [ndwi,f_name_out] = atm_CAMBOT_calc_NDWI_L1B(f_name_cambot,save_geotiff);\n',...
        '\tSYNTAX: [ndwi]            = atm_CAMBOT_calc_NDWI_L1B(f_name_cambot);\n'])
end

% check if input file exists
if (ischar(f_name_cambot) && exist(f_name_cambot) == 0)
        warndlg({'Input file not found!';'Script aborted.'},'!! Warning !!');
        error('atm_CAMBOT_calc_NDWI_L1B:file_chk', '\n\tInput file not found. Script aborted.')
end

% determine if geotiff output is desired
if (length(varargin) >=1 )
    save_geotiff = varargin{1};
    if (save_geotiff ~= 0 && save_geotiff ~= 1)
        error('atm_CAMBOT_calc_NDWI_L1B:save_geotiff', '\n\tERROR: save_geotiff must be either 0 or 1.')
    end
elseif(nargin == 1 || exist('save_geotiff') == 0) % save_geotiff has not been set due to nargin == 1
    save_geotiff = 0;
end


%% import CAMBOT natural color geotiff image

info_tiff       = geotiffinfo(f_name_cambot);    % read parameters needed for exporting geolocated geotiff file
[cambot_rgb,R0] = readgeoraster(f_name_cambot);  % read raster and map cell reference RO for exporting geotiff

%% calculate NDWI_ice = (blue - red)/(blue + red); 
%  NDWI_ice values < 0 sometimes indicate clouds, rain, and snow 

cambot_rgb = single(cambot_rgb);  % needs to be converted to single precision for NDWI_ice calculation
ndwi = (cambot_rgb(:,:,3) - cambot_rgb(:,:,1))./(cambot_rgb(:,:,3) + cambot_rgb(:,:,1));


%% export geotiff if desired

if (save_geotiff == 1)
    
    % build output file name
    [f_path,f_name,f_ext] = fileparts(f_name_cambot);
    f_name_out = [f_path filesep f_name '_NDWI' f_ext];
    
    % write NDWI_ice geotiff file using LZW compression  
    tags.Compression = Tiff.Compression.LZW;
    geotiffwrite(f_name_out,ndwi,R0,'GeoKeyDirectoryTag',info_tiff.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags',tags);
    
    listing = dir(f_name_out);
    fprintf('Exported %4.1f Mbytes in %s\n',listing.bytes/1024/1024,f_name_out);
    
end

%% prepare output parameters

n_outputs = nargout;
varargout = cell(1,n_outputs);
varargout{1} = ndwi;

if n_outputs == 2
    varargout{2} = f_name_out;
end

end