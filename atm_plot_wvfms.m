function [fig_wvfm] = atm_plot_wvfms(atm_wvfm,varargin)
%
% SUMMARY:           ATM_PLOT_WVFMS plots ATM waveform data stored in a stuctured array created by ATM_WVFM_READER
%                   
%                    The number of input arguments must be between 2 and 5 with the following options:
% SYNTAX:            plt_wvfms(wvfm_struct,shot_list);
% SYNTAX:            plt_wvfms(wvfm_struct,shot_list,t_min);
% SYNTAX:            plt_wvfms(wvfm_struct,shot_list,t_min,t_max);
% SYNTAX:            plt_wvfms(wvfm_struct,shot_list,t_min,t_max,back_ground_color);
% 
% wvfm_struct:       First input argument must be a stuctured array created by ATM_WVFM_READER
%
% shot_list:         Second input argument must be a row or column vector of shot indices to be plotted.
%
% t_min, t_max:      Laser trigger time in nano seconds. PLT_WVFMS plots only range gates with t >= t_min and t <= t_max
%                    t_max can be let unset or set to Inf to include all range gates after t_min
%                    The two parameters can be used to select window reflections, transmit pulses or returns from targets
%                    The time of the window reflection and transmit pulse changes with installations and needs to be
%                    adjusted accordinly.
%                    A good starting point to only plot windwow reflections is
%                        t_min = 0, t_max = 80
%                    Plot only transmit pulses:
%                        t_min = 80, t_max = 200
%                    Plot all target returns without window reflections or transmit pulses
%                        t_min = 200, tmax = Inf (or left blank)
%
%                    The time windows for window reflections, transmit pulse and return pulsea to use for the particular data file
%                    are campaign and instrument dependent and are listed in the HDF5 fields /waveforms/twv/ancillary_data/rx_start, 
%                    ./tx_start and ./tx_end in number of digitizer samples.
%
% back_ground_color: can be 'w' (white) or 'b' (black). If not set PLT_WVFMS create a figure with black background by default. 
% 
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Author:            Michael Studinger, NASA Goddard Space Flight Center, Greenbelt MD, USA.
% Version:           3.04 - March 11, 2019
% See also:          ATM_WVFM_READER
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------


%% verify input arguments

% first check number of input arguments
if ~any(nargin == [2 3 4 5])
    error('plt_wvfms:t_min_t_max_plot', ['\n\tERROR: The number of input arguments must be between 2 and 5 with the following options:\n\n', ...
        '\tSYNTAX: plt_wvfms(wvfm_struct,shot_list);\n',...
        '\tSYNTAX: plt_wvfms(wvfm_struct,shot_list,t_min);\n',...
        '\tSYNTAX: plt_wvfms(wvfm_struct,shot_list,t_min,t_max);\n',...
        '\tSYNTAX: plt_wvfms(wvfm_struct,shot_list,t_min,t_max,back_ground_color);\n'])
end

% check if first input argument is a MATLAB structured array and was likely created from a function call to atm_wfvm_reader
if ~(isstruct(atm_wvfm)) % check if structured array
    error('plt_wvfms:struct1','\n\tERROR: First input argument must be a stuctured array created by atm_wvfm_reader.\n');
elseif ~(isfield(atm_wvfm,'shot_id')) % check if atm_wvfm_reader has created the struct
    error('plt_wvfms:struct2','\n\tERROR: First input argument appears to be a stuctured array but was not created by atm_wvfm_reader.\n');
end

% check if second input argument is a valid shot list
if (nargin >= 2 && isnumeric(varargin{1}) &&  isvector(varargin{1})) % check if this is a numeric array
    shot_list = varargin{1};
    if (isrow(shot_list) == 1)          % determine the number of shots for a row vector
        n_shots = size(shot_list,2);
    elseif (iscolumn(shot_list) == 1)   % determine the number of shots for a column vector
        shot_list = shot_list';
        n_shots = size(shot_list,2);
    end
else
    error('plt_wvfms:shot_list1','\n\tERROR: Second input argument must be a row or colum vector of shot indices to be plotted.\n');
end

% create warning if shot list contains large number of shots, but continue plotting
if (n_shots >= 200)
    warning('plt_wvfms:n_shots', '\n\tWARNING: This function is intended to plot a small number of waveforms.\n')
end

% now check for t_min and t_max
if (nargin == 3 && isnumeric(varargin{2})) % only t_min set
    t_min = varargin{2};
    t_max = Inf; % Inf returns the IEEE® arithmetic representation for positive infinity if t_max is not set by input argument
elseif (nargin >= 4 && isnumeric(varargin{2}) && isnumeric(varargin{3})) % this is for t_min and t_max
    t_min = varargin{2};
    t_max = varargin{3};
else
    error('plt_wvfms:t_min_t_max_plot', ['\n\tERROR: The number of input arguments must be between 2 and 5 with the following options:\n\n', ...
        '\tSYNTAX: plt_wvfms(wvfm_struct,shot_list);\n',...
        '\tSYNTAX: plt_wvfms(wvfm_struct,shot_list,t_min);\n',...
        '\tSYNTAX: plt_wvfms(wvfm_struct,shot_list,t_min,t_max);\n',...
        '\tSYNTAX: plt_wvfms(wvfm_struct,shot_list,t_min,t_max,back_ground_color);\n'])
end
    
% now check for background color
if (nargin == 5 && ischar(varargin{4})) % determine background color if set
    back_ground_color = varargin{4};
else % default setting if not set
    back_ground_color = 'b';
end

% check for valid background color setting
if ~(strcmp(back_ground_color,'w') || strcmp(back_ground_color,'b'))
    warning('plt_wvfms:color1', ['\n\tWARNING: background color should be b or w or not set.',...
        '\n\tUsing default background color (black).'])
    back_ground_color = 'b';
end

% verify t_min < t_max
if (t_min >= t_max)
    error('plt_wvfm:t_start_t_end', '\n\tERROR: t_min must be < t_max.')
end

%% plot 
edge = 50; scrsz = get(0,'ScreenSize');

fig_wvfm = figure('Position',[edge edge (scrsz(3)-4*edge) (scrsz(4)-3*edge)]); 

if strcmp(back_ground_color,'b') % set figure background to black 
    colordef(fig_wvfm,'black');
end

for n = 1:n_shots
    shot_nr = shot_list(n);
    for i = 1:atm_wvfm.n_gates(shot_nr)
        if( atm_wvfm.shots(shot_nr).wf(i).t(1) >= t_min && atm_wvfm.shots(shot_nr).wf(i).t(1) <= t_max)
            shot_id_label = sprintf('%d',atm_wvfm.shot_id(shot_nr));
            plot(atm_wvfm.shots(shot_nr).wf(i).t,atm_wvfm.shots(shot_nr).wf(i).w,'-','DisplayName',shot_id_label,'LineWidth',1.75); hold on;
            clear shot_id_label;
        end
    end
end

h = legend('show','Location','northeastoutside');
title(h,'Unique Shot Identifier');

grid on;
xlabel('Laser Trigger Time [ns]');
ylabel('Signal Amplitude [counts]');

if (atm_wvfm.info.ground_test_data == 0)
    
    if atm_wvfm.info.laser_wavelength_nm == 532
        title(sprintf('%s %d nm: %d kHz, %3.1f ns, %.1f°',atm_wvfm.info.sensor,atm_wvfm.info.laser_wavelength_nm,...
            atm_wvfm.info.laser_prf_hz/1000,atm_wvfm.info.pulse_width_fwhm_ns,atm_wvfm.info.laser_off_nadir_angle),'Color','g');
    elseif atm_wvfm.info.laser_wavelength_nm == 1064
        title(sprintf('%s %d nm: %d kHz, %3.1f ns, %.1f°',atm_wvfm.info.sensor,atm_wvfm.info.laser_wavelength_nm,...
            atm_wvfm.info.laser_prf_hz/1000,atm_wvfm.info.pulse_width_fwhm_ns,atm_wvfm.info.laser_off_nadir_angle),'Color','r');
    else
        title(sprintf('%s %d nm: %d kHz, %3.1f ns, %.1f°',atm_wvfm.info.sensor,atm_wvfm.info.laser_wavelength_nm,...
            atm_wvfm.info.laser_prf_hz/1000,atm_wvfm.info.pulse_width_fwhm_ns,atm_wvfm.info.laser_off_nadir_angle),'Color','w');
    end
    
end
set(gca,'FontSize',16);


end

