function [atm_wvfm] = atm_wvfm_ground_test_reader(f_name_inp,varargin)
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% NOTE:             TEMPORAL and SPATIAL searches are disabled for reading ground test waveform data.
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
%
% SUMMARY:          The ATM_WVFM_GROUND_TEST_READER function reads ATM waveform data from an ATM HDF5 ground test waveform files and imports the data into a 
%                   structured array (struct) in MATLAB. This function is a version of the atm_wvfm_reader function with spatial and temporal search capabilities disabled
%                   because ground test waveform data do not contain spatial or temporal information.
%
% SYNTAX:           atm_wvfm = atm_wvfm_ground_test_reader(f_name_inp);                                              % silent conversion of entire file   
%                   atm_wvfm = atm_wvfm_ground_test_reader(f_name_inp, verbose);                                     % conversion of entire file with optional console output
%
% EXAMPLE:          Example of a valid function call with command window/console output:
%                   atm_wvfm = atm_wvfm_ground_test_reader('ILATM1B_20170510_132857.atm6AT6.h5',1);
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% INPUT PARAMETERS
% 
%   f_name_inp      File name of HDF5 ATM ground test waveform file to be read including pathname.
% 
%   verbose         Must be 0 or 1. 1 displays some basic parameters on the console. 0 is silent.
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% OUTPUT PARAMETERS
%
%   atm_wvfm        Varible name of MATLAB structured array (struct) containing waveform data.
%                   Metadata and processing information such as search criteria, instrument configuration and data format version as included as well.
%                   atm_wvfm must be a valid MATLAB variable name.  
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Author:           Michael Studinger, NASA Goddard Space Flight Center, Greenbelt MD, USA.
% Version:          4.01 - February 24, 2022
% See also:         atm_wvfm_info.m, atm_wvfm_reader.m
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------

%% parse input parameters and do basic error checking

% first determine if this is a groundtest file

grd_test = h5read(f_name_inp,'/aircraft/latitude');
if (length(grd_test) == 0)
    %fprintf('\nWarning: This appears to be a ground test file. No spatial or temporal search is supported.\n');
    opts = struct('WindowStyle','modal','Interpreter','tex');
    h = warndlg({'\color{red}This appears to be a ground test file. No spatial or temporal search is supported.';'Press OK to continue.'},'!! Warning !!', opts);
    uiwait(h,1); % automatically close warning dialog after 5 seconds and continue
    %close(h);
    fprintf('Importing ground test data...\n');
    GROUND = 1;
else
    GROUND = 0;
end
clear grd_test;

% check if input file exists
if (nargin >= 1 && ischar(f_name_inp))
    if (exist(f_name_inp) == 0)
        warndlg({'Input file not found!';'Script aborted.'},'!! Warning !!');
        error('atm_wvfm_reader:file_chk', '\n\tInput file not found. Script aborted.')
    end
end

% determine if console output is desired
if (length(varargin) >= 1)
    verbose = varargin{1};
    if (verbose ~= 0 && verbose ~= 1)
        error('atm_wvfm_info:verbose', '\n\tERROR: verbose must be either 0 or 1.')
    end
elseif(nargin == 1 || exist('verbose') == 0) % verbose has not been set due to nargin == 1
    verbose = 0;
end

%% determine version of ATM waveform data format

info = h5info(f_name_inp);
if (size(info.Groups,1) == 7)
    if (isfield(info,'Groups') && strcmp(info.Groups(7).Groups(1).Name,'/waveforms/twv'))
        data_fmt = 2; % current ATM waveform format stored in subgroup /waveforms/twv
        data_fmt_str = 'ATM Waveform Format Vers. 2.0';
    elseif (isfield(info,'Groups') && strcmp(info.Groups(8).Groups.Name,'/waveforms/vld'))
        data_fmt = 1; % legacy ATM waveform format stored in subgroup /waveforms/vld
        data_fmt_str = 'ATM Waveform Format Vers. 1.0';
        error('atm_wvfm_reader:wvfm_fmt_vers', '\n\tERROR: ATM waveform format version 1 not yet implemented.')
    else
        error('atm_wvfm_reader:wvfm_fmt_vers', '\n\tERROR: Cannot determine ATM waveform format version from HDF5 file.')
    end
    
elseif (size(info.Groups,1) == 6) % IR file withouth /footprint directory
    if (isfield(info,'Groups') && strcmp(info.Groups(6).Groups(1).Name,'/waveforms/twv'))
        data_fmt = 0; % current ATM waveform format stored in subgroup /waveforms/twv
        data_fmt_str = 'ATM Waveform Format Vers. 2.0 (IR)';
    end
else
    error('atm_wvfm_info:input_struct','\n\tERROR: ATM waveform struct must have 7 groups.\n')
end

clear info;

%% read data and set up MATLAB serial time tags for individual laser shots

% read time tags (seconds of day past midnight) of laser trigger time
shots_seconds_of_day = h5read(f_name_inp,'/waveforms/twv/shot/seconds_of_day');


%% determine the record/shot and range gate information necessary to reassemble the waveforms
% information is stored in different subgroups depending on data format type: 1 (ATM legacy) or 2 (current waveform)
t1 = tic;

if (data_fmt == 2 || data_fmt == 0)
        
    % sampling interval
    sample_int_ns = h5read(f_name_inp,'/waveforms/twv/ancillary_data/sample_interval'); % sampling interval in nano seconds
    
    % laser shots
    shot_gate_count = h5read(f_name_inp,'/waveforms/twv/shot/gate_count');   % number of gates in shot
    shot_gate_start = h5read(f_name_inp,'/waveforms/twv/shot/gate_start');   % 1-based record number of waveform's first gate in shot (1, 4, 7,...)
    shot_identifier = h5read(f_name_inp,'/waveforms/twv/shot/number');       % unique shot identifier
    shot_gate_xmt   = h5read(f_name_inp,'/laser/gate_xmt');                  % number of 1-based range gate in shot/record with transmit waveform (typically 1 or 2)
    shot_gate_rcv   = h5read(f_name_inp,'/laser/gate_rcv');                  % number of 1-based range gate in shot/record with return waveform (typically 2, 3 or higher)
    cal_range       = h5read(f_name_inp,'/laser/calrng');                    % calibrated slant range from ATM origin to surface
    
    rx_width = h5read(f_name_inp,'/laser/pulse_width');     
    rx_sigstr = h5read(f_name_inp,'/laser/rcv_sigstr'); 
    
    % gate information ~ satu count not yet used. 
    gate_wvfm_start  = h5read(f_name_inp,'/waveforms/twv/gate/wvfm_start');  % int32 pointer to 1-based first sample in waveform gate (1, 97, 193, ...)
    gate_wvfm_length = h5read(f_name_inp,'/waveforms/twv/gate/wvfm_length'); % int32
    gate_position    = h5read(f_name_inp,'/waveforms/twv/gate/position');    % int32 number of samples after laser trigger (n = 0 marks trigger)
    gate_satu_cnt    = h5read(f_name_inp,'/waveforms/twv/gate/pulse/sat_count');  % int32 number of samples in gate that are saturated (255 counts for 8 bit digitizer)
    
    % added June 27, 2018 & Dec 18, 2018
    gate_pulse_width = h5read(f_name_inp,'/waveforms/twv/gate/pulse/width');    % width of pulse (number of samples) based on a threshold of 35% of peak
    gate_n_pulses    = h5read(f_name_inp,'/waveforms/twv/gate/pulse/count');    % number of pulses in gate (= number of threshold crossings divided by two)
    gate_area        = h5read(f_name_inp,'/waveforms/twv/gate/pulse/area');    % area of waveform pulse above noise floor (in counts * nanoseconds)
    
    % some useful parameters, added July 2, 2018
    gate_xmt         = h5read(f_name_inp,'/laser/gate_xmt');                    % gate number of transmit pulse
    gate_rcv         = h5read(f_name_inp,'/laser/gate_rcv');                    % gate number of primary receive pulse
    
    % range bins for faster import
    wvfm_amplitude = h5read(f_name_inp,'/waveforms/twv/wvfm/amplitude');     % uint8

    
elseif (data_fmt == 1)
    error('atm_wvfm_reader:data_fmt', '\n\tERROR: vld data format not yet supported.')
end
    

%% extract range gates and arrange them into a MATLAB structure array

% populate the waveform struct - this is for some reason not much faster than doing it inside the loop

shot_list_search = 1:length(shot_identifier);

atm_wvfm = struct('shot_id',shot_identifier(shot_list_search),...
    'n_gates',shot_gate_count(shot_list_search),...
    'shot_gate_start', shot_gate_start(shot_list_search),...
    'shot_gate_xmt',shot_gate_xmt(shot_list_search),...
    'shot_gate_rcv',shot_gate_rcv(shot_list_search),...
    'cal_range',cal_range(shot_list_search),...
    'rx_width',rx_width(shot_list_search),...
    'rx_sigstr',rx_sigstr(shot_list_search),...
    'gate_xmt',gate_xmt(shot_list_search),...
    'gate_rcv',gate_rcv(shot_list_search));   
    
tStart = tic;

for i = 1:length(shot_list_search)
    
    record_nr = shot_list_search(i);
    
    n_gates = shot_gate_count(record_nr); % determine the number or range gates
        
    for k = 1:n_gates
        
        % determine pointers/indices of desired shot/record number and accompanying range gates
        indx = shot_gate_start(record_nr) + k - 1;
        atm_wvfm.shots(i).wf(k).w = wvfm_amplitude(gate_wvfm_start(indx):gate_wvfm_start(indx) + gate_wvfm_length(indx) - 1); 

        % populate array of laser trigger time tags in units of nano seconds - convenient, but results in large data volume, when used for entire files
        tw_tmp = 0:sample_int_ns:double(gate_wvfm_length(indx))-1;
        tw_tmp1 = tw_tmp(1:gate_wvfm_length(indx));
        atm_wvfm.shots(i).wf(k).t = double(gate_position(indx))*sample_int_ns + tw_tmp1;
        
        % add number of saturated samples and pulse width
        atm_wvfm.shots(i).wf(k).n_sat = gate_satu_cnt(indx);
        atm_wvfm.shots(i).wf(k).pulse_width = gate_pulse_width(indx);
        atm_wvfm.shots(i).wf(k).gate_n_pulses = gate_n_pulses(indx); 
        atm_wvfm.shots(i).wf(k).area = gate_area(indx);
        
    end
    
    clear record_nr indx tw_*;
    
end

t_elapsed = toc(tStart);

%% create info field containing information about the subsetting/HDF5 ingest and basic metadata about the HDF5 input file

atm_wvfm.info.f_name                 = char(f_name_inp);
atm_wvfm.info.data_fmt_str           = data_fmt_str;
atm_wvfm.info.sampling_int_ns        = sample_int_ns;
atm_wvfm.info.laser_prf_hz           = h5read(f_name_inp,'/ancillary_data/instrument/laser_prf');
atm_wvfm.info.pulse_width_fwhm_ns    = h5read(f_name_inp,'/ancillary_data/instrument/laser_pulse_width');
atm_wvfm.info.laser_wavelength_nm    = h5read(f_name_inp,'/ancillary_data/instrument/laser_wavelength');
atm_wvfm.info.rx_start_ns            = h5read(f_name_inp,'/waveforms/twv/ancillary_data/rx_start')*sample_int_ns;
atm_wvfm.info.tx_start_ns            = h5read(f_name_inp,'/waveforms/twv/ancillary_data/tx_start')*sample_int_ns;
atm_wvfm.info.tx_end_ns              = h5read(f_name_inp,'/waveforms/twv/ancillary_data/tx_end')*sample_int_ns;
atm_wvfm.info.rx_start_samples       = h5read(f_name_inp,'/waveforms/twv/ancillary_data/rx_start'); % 1-based starting bound for transmit waveform search in digitizer samples
atm_wvfm.info.tx_start_samples       = h5read(f_name_inp,'/waveforms/twv/ancillary_data/tx_start'); % 1-based starting bound for receive waveform search in digitizer samples
atm_wvfm.info.tx_end_samples         = h5read(f_name_inp,'/waveforms/twv/ancillary_data/tx_end');   % 1-based endpoint for transmit waveform search in digitizer samples
atm_wvfm.info.date_processed         = datestr(now);
atm_wvfm.info.m_file_used            = mfilename('fullpath');
atm_wvfm.info.user                   = getenv('UserName');
atm_wvfm.info.atm_processing_version = h5read(f_name_inp,'/ancillary_data/documentation/version');
atm_wvfm.info.atm_header_text        = h5read(f_name_inp,'/ancillary_data/documentation/header_text');
atm_wvfm.info.input_parameters       = varargin;
atm_wvfm.info.ground_test_data       = GROUND;

%% display processing information in MATLAB Command Window or Console if desired

if (verbose == 1)
    
    [~,name,ext] = fileparts(f_name_inp); 
    fprintf('\n-----------------------------------------------------------------------\n');
    fprintf('Imported file %s\n',[name ext]);
    fprintf('-----------------------------------------------------------------------\n');
    fprintf('Processing time (sec):                    %6.2f\n',t_elapsed);
    fprintf('Number of laser shots in ATM HDF5 file: %8d\n',length(cal_range));
    fprintf('-----------------------------------------------------------------------\n');
    
end

end
