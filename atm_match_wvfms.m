function [indx_gr, indx_ir, varargout] = atm_match_wvfms(f_name_gr,varargin)
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% SUMMARY:          The ATM_MATCH_WVFMS function determines matching laser shots from the ATM dual-color (532 and 1064 nm) narrow transceiver (ATM-T7).
%                   Both NIR (1064 nm) and green (532 nm) wavelength measurements are acquired by a single laser and transceiver, but by separate data systems.
%                   The measurements are co-located spatially and synchronized temporally. The NIR and green laser shots are stored in separate HDF5 data files.
%
%                   The green data product is ILNSAW1B (https://nsidc.org/data/ILNSAW1B) and the IR data product is ILNIRW1B (https://nsidc.org/data/ilnirw1b).
%                   The data records from the two data systems need to be matched, as occasionally one or the other system will not exceed the trigger threshold
%                   and therefore won't record data for that laser shot. 
%                   This co-registration can be performed with the ATM_MATCH_WVFMS function for finding matching shots in both data products.
%
% DEPENDENCY:       ATM_MATCH_WVFMS requires the binary search function from Avi Ziskind written in C and compiled as a binary MEX file.
%                   The source code can be downloaded from MATLAB's File Exchange page:
%                   https://www.mathworks.com/matlabcentral/fileexchange/30484-fast-binary-search
%                   To compile the C code run: mex binarySearch.c
%                   The user can determine the platform specific file-name extension of the binary MEX file with the MATLAB command "mexext".
%                   On Windows 64 bit systems the MEX file is binarySearch.mexw64. The MEX file name on Mac Platforms is binarySearch.mexmaci64.
%                   For 64 bit Linux platforms the compiled MEX file is binarySearch.mexa64.
%
% SYNTAX:           [indx_gr, indx_ir] = atm_match_wvfms(f_name_inp);                             % silent co-registration of NIR and green data files   
%                   [indx_gr, indx_ir] = atm_match_wvfms(f_name_inp, verbose);                    % co-registration with optional console output
%                   [indx_gr, indx_ir] = atm_match_wvfms(f_name_inp, verbose, path_ir_file);      % pathname to NIR files if different than green files
%
% EXAMPLES:         Example of a valid function call for matching NIR and green files and command window/console output:
%                   [indx_gr, indx_ir] = atm_match_wvfms('ILNSAW1B_20181104_164800.atm6DT7.h5',1);
%                   [indx_gr, indx_ir] = atm_match_wvfms_v1e('D:\atm\green_ir_matching\ILNSAW1B_20181104_164800.atm6DT7.h5',1);
%                   [indx_gr, indx_ir, info] = atm_match_wvfms_v1e('D:\atm\green_ir_matching\ILNSAW1B_20181104_164800.atm6DT7.h5',1,'D:\atm\green_ir_matching\NIR');
%                   [indx_gr, indx_ir, info] = atm_match_wvfms_v1e('D:\atm\green_ir_matching\ILNSAW1B_20181104_164800.atm6DT7.h5',1);
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% INPUT PARAMETERS
% 
% f_name_inp        File name of HDF5 ATM ILNSAW1B green waveform filr to be used for co-registration
% verbose           Must be 0 or 1. 1 displays some basic parameters on the console. 0 is silent.
% path_ir_file      pathname to NIR file if the corresponding NIR data file (ILNIRW1B) is in a different folder than the green file (ILNSAW1B)  
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% OUTPUT PARAMETERS
%
% indx_gr, indx_ir  indices for green and NIR HDF5 files that result in co-registered, matching NIR and green laser shots
%                   indx_gr and indx_ir must be valid MATLAB variable names.
% info              [optional] MATLAB struct with matching IR file name, metadata and statistics created during processing
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Author:           Michael Studinger, NASA Goddard Space Flight Center, Greenbelt MD, USA.
% Version:          1.1 - March 3, 2022
% Compatibility:    tested with MATLAB Version 9.11.0.1769968 (R2021b) on Windows 10 Pro 64 bit (Version 21H2).
% See also:         atm_wvfm_reader.m, atm_wvfm_info.m and https://www.mathworks.com/matlabcentral/fileexchange/30484-fast-binary-search
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------

%% parse input parameters and do basic error checking

% check number of input arguments
if ~any(nargin == [1 2 3])
    error('atm_match_wvfms:nargin', ['\n\tERROR: The number of input arguments must be 1, 2, or 3:\n\n', ...
        '\tSYNTAX: [indx_gr, indx_ir] = atm_match_wvfms(f_name_gr);\n',...
        '\tSYNTAX: [indx_gr, indx_ir] = atm_match_wvfms(f_name_gr, verbose);\n',...
        '\tSYNTAX: [indx_gr, indx_ir] = atm_match_wvfms(f_name_gr, verbose, path_ir_file);\n'])
end

% check number of output arguments - alternatively use nargoutchk(2,2)
if ~any(nargout == [2 3])
    error('atm_match_wvfms:nargout', ['\n\tERROR:  Number of output arguments must be 2 or 3:\n', ...
        '\tSYNTAX: [indx_gr, indx_ir] = atm_match_wvfms(f_name_gr);\n',...
        '\tSYNTAX: [indx_gr, indx_ir] = atm_match_wvfms(f_name_gr, verbose);\n',...
        '\tSYNTAX: [indx_gr, indx_ir] = atm_match_wvfms(f_name_gr, verbose, path_ir_file);\n'])
end

% determine if console output is desired
if (length(varargin) >= 1)
    verbose = varargin{1};
    if (verbose ~= 0 && verbose ~= 1)
        error('atm_match_wvfms:verbose', '\n\tERROR: verbose must be either 0 or 1.')
    end
elseif(nargin == 1 || exist('verbose') == 0) % verbose has not been set due to nargin == 1
    verbose = 0;
end

% check if input file exists
if (nargin >= 1 && ischar(f_name_gr))
    if (exist(f_name_gr) == 0)
        disp(f_name_gr)
        warndlg({'Input file not found!';'Script aborted.'},'!! Warning !!');
        error('atm_match_wvfms:file_chk', '\n\tInput file with green shots not found. Script aborted.')
    end
end

%% build file names. make search for IR independent from future data system changes and transceiver versions, which will result in different file names.

if (length(varargin) == 2 & exist(varargin{2}) == 7) % make sure path name exists and is a folder
    f_path_ir = char(varargin{2});
end

f_name_ir = make_ir_f_name(f_name_gr);

%% set up struct for loops

f_names(1).name = f_name_gr;
f_names(2).name = f_name_ir;
f_names(1).laser = 'gr';
f_names(2).laser = 'ir';

%% check if both input files exist

for i = 1:length(f_names)
    
    f_info = dir(f_names(i).name);
    
    if (exist(f_names(i).name) == 0)
        warndlg({'Input file not found!';'Script aborted.'},'!! Warning !!');
        error('atm_match_wvfms:file_chk', '\n\tInput file %s not found. Script aborted.',f_names(i).name);
    elseif (exist(f_names(i).name) == 2 & verbose == 1)
        fprintf('%s checked: %05.1f MB\n',f_names(i).name,f_info.bytes/1024/1024);
    end
    
    clear f_info;
end

%% check if files appear to be a valid ATM HDF5 files from NSIDC

for i = 1:length(f_names)
    
    info = h5info(f_names(i).name);
    
    if (size(info.Groups,1) == 7 & f_names(i).laser == 'gr') % green file
        if (isfield(info,'Groups') && strcmp(info.Groups(7).Groups(1).Name,'/waveforms/twv'))
            data_fmt = 2; % current ATM waveform format stored in subgroup /waveforms/twv
            data_fmt_str = 'ATM Waveform Format Vers. 2.0';
            
            laser_w_length_gr = h5read(f_name_gr,'/ancillary_data/instrument/laser_wavelength');
            if laser_w_length_gr ~= 532
               error('atm_wvfm_match:wvfm_fmt_vers', '\n\tERROR: Laser wavelength is not 532 nm. Script aborted.') 
            end
        elseif (isfield(info,'Groups') && strcmp(info.Groups(8).Groups.Name,'/waveforms/vld'))
            data_fmt = 1; % legacy ATM waveform format stored in subgroup /waveforms/vld
            data_fmt_str = 'ATM Waveform Format Vers. 1.0';
            error('atm_wvfm_match:wvfm_fmt_vers', '\n\tERROR: ATM waveform format version 1 not yet implemented.')
        else
            error('atm_wvfm_match:wvfm_fmt_vers', '\n\tERROR: Cannot determine ATM waveform format version from HDF5 file.')
        end
        
    elseif (size(info.Groups,1) == 6 & f_names(i).laser == 'ir') % IR file withouth /footprint directory
        if (isfield(info,'Groups') && strcmp(info.Groups(6).Groups(1).Name,'/waveforms/twv'))
            data_fmt = 0;  % the current ATM waveforms are stored in subgroup /waveforms/twv
            data_fmt_str = 'ATM Waveform Format Vers. 2.0 (IR)';
            
            laser_w_length_ir = h5read(f_name_ir,'/ancillary_data/instrument/laser_wavelength');
            if laser_w_length_ir ~= 1064
               error('atm_wvfm_match:wvfm_fmt_vers', '\n\tERROR: Laser wavelength is not 532 nm. Script aborted.') 
            end
        end
    else
        error('atm_wvfm_match:input_struct','\n\tERROR: ATM IR waveform input file must have 6 groups.\n')
    end
    
    clear info;
    
end

%% determine tolerance threshold based on sampling interval

% search tolerance for finding matching shots (should be less than 50 usec, i.e., half the laser interval)
% tolerance = 40e-6; % 40 us (most matches are within 10 usec)

sample_int_ns = h5read(f_name_gr,'/waveforms/twv/ancillary_data/sample_interval'); % sampling interval in nano seconds
sample_int = sample_int_ns*1E-9;      % sampling interval in seconds
tolerance = 100000*2*sample_int/1.25; % not clear how this relates to sampling interval

%% read utc time stamps and make sure they are sorted. this will also identify identical time tags, which should not be in the data.

utc_gr = h5read(f_name_gr,'/time/seconds_of_day');    
utc_ir = h5read(f_name_ir,'/time/seconds_of_day');

% make sure timetags are sorted - this should be the case
if (issorted(utc_gr) == 0)
    error('atm_wvfm_match:utc_gr_sorted', '\n\tERROR: UTC time tags (532 nm) are not sorted. Script aborted.')
elseif (issorted(utc_ir) == 0)
    error('atm_wvfm_match:utc_ir_sorted', '\n\tERROR: UTC time tags (1064 nm) are not sorted. Script aborted.')
end

%% search for matching laser shots

tStart = tic;

% allocate arrays for faster execution
indx_ir = (1:length(utc_ir))';
indx_gr = zeros(length(utc_ir),1); indx_gr = indx_gr./0; % set all elements to NaN 
delta_t = zeros(length(utc_ir),1); delta_t = delta_t./0; % set all elements to NaN

% this loop is the meat of all of it
for i=1:length(utc_ir)
    indx_gr(i) = binarySearch(utc_gr, utc_ir(i),[],2); % returns the position closest to utc_ir(i) in utc_gr 
end

% vectorize delta_t computation and determine which IR laser shots have no corresponding green shots
delta_t = utc_gr(indx_gr) - utc_ir;
no_match = find(abs(delta_t) > tolerance);

% set up indices for both green and IR
indx_match = find(abs(delta_t) < tolerance);
indx_gr = indx_gr(indx_match);
indx_ir = indx_ir(indx_match);

% caluclate some statistics
n_matching_shots = length(utc_ir)-length(no_match);
n_nomatch_shots = length(no_match);
delta_t_matching_sec = mean(delta_t(indx_match));
std_t_matching_sec   = std(delta_t(indx_match));

t_elapsed = toc(tStart);


%% populate info struct with processing information if desired

if (nargout == 3)
    
   info.orig_shots_gr = length(utc_gr); 
   info.orig_shots_ir = length(utc_ir);
   info.n_matching_shots = n_matching_shots;
   info.n_nomatch_shots  = n_nomatch_shots;
   info.search_tolerance_sec = tolerance;
   info.delta_t_matching_sec = delta_t_matching_sec;
   info.std_t_matching_sec = std_t_matching_sec;
   info.f_name_gr = f_name_gr;
   info.f_name_ir = f_name_ir;
   
   varargout{1} = info;
   
end

%% display processing information in MATLAB Command Window or Console if desired

if (verbose == 1)

    fprintf('\n-----------------------------------------------------------------------\n');
    fprintf('Search for matchin IR and green laser shots.  Tolerance = %.7f sec\n',tolerance);
    fprintf('-----------------------------------------------------------------------\n');
    fprintf('Processing time (sec)           : %6.2f\n',t_elapsed);
    fprintf('Number of laser shots in GR file: %8d\n',length(utc_gr));
    fprintf('Number of laser shots in IR file: %8d\n',length(utc_ir));
    fprintf('Number of matching laser shots  : %8d\n',n_matching_shots);
    fprintf('Green shots without matching IR :   %d (%.1f%%)\n',n_nomatch_shots,100*n_nomatch_shots/(length(utc_ir)-n_nomatch_shots));
    fprintf('Mean time difference [secs]     :   %9.7f\n',delta_t_matching_sec);
    fprintf('Std time difference  [secs]     :   %9.7f\n',std_t_matching_sec);

    fprintf('-----------------------------------------------------------------------\n');
    
end

%% function to build corresponding IR file name based on green file name independent from future data system changes and transceiver versions

function f_ir = make_ir_f_name(f_gr)

    % get path name, file name, and extension for the specified file
    [f_path, f_name, f_ext] = fileparts(f_gr);
    % for NIR files stored in a different folder than the green files
    if (length(varargin) == 2 & exist(varargin{2}) == 7) % make sure path is an existing directory
        f_path = f_path_ir;
    end       
    % first, replace the green data product designation ILNSAW1B with IR data product designation ILNIRW1B
    f_name = strrep(f_name,'ILNSAW1B','ILNIRW1B');
    % remove .atm6DT7 or any other transeiver and data system version at the end
    [tmp, f_name_tmp, tmp] = fileparts(f_name); clear tmp;
    % search for IR file independent of future changes in data system and transceiver version
    listing = dir(sprintf('%s%s%s.*.h5',f_path, filesep, f_name_tmp));  clear tmp f_name_tmp f_path f_ext;

    if (length(listing) == 1) % found corresponding IR file, make sure it is one unique file
        f_ir = [listing.folder filesep listing.name];
    elseif (length(listing) ~= 1)
        error('atm_match_wvfms:file_chk', '\n\tUnable to locate unique IR file for %s. Script aborted.',f_gr);
    end

    clear listing;
    
end % function f_ir = make_ir_f_name(f_gr)


end % function [indx_gr, indx_ir] = atm_match_wvfms(f_name_gr,varargin)
