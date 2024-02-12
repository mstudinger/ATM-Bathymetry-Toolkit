function [w_ti, t_i, t_c] = atm_centroid_tracker(t,w,threshold)
%
% SUMMARY:            The ATM_CENTROID_TRACKER function is a MATLAB implementation of ATM's way to estimate two-way travel times using the centroid of the 
%                     transmit and return waveforms. The centroid tracker is using the part of the waveform that is above 15% of the maximum amplitude above
%                     the noise floor. The noise floor is estimated using the median of the first 21 samples of the waveform.
%                     If used with a threshold of 0.5 the function can be used to estimate the full width at half maximum (FWHM) of
%                     the pulse. The output is the waveform for the part above the threshold (t_i, w_i) and the time estimate of the centroid (t_c).
%
% SYNTAX:             [w_ti, t_i, t_c] = atm_centroid_tracker(t,w,threshold);
%                 
% --------------------------------------------------------------------------------------------------------------------------------------------------------------
%
% INPUT PARAMETERS
% 
%   t                 Waveform time tags in nanoseconds as a column vector. When atm_wvfm_reader.m is used for importing the ATM waveform data 
%                     t = [atm.shots(k).wf(i).t];
%
%   w                 Corresponding waveform amplitudes for t as a column vector. w and t must have the same length. When atm_wvfm_reader.m is used for 
%                     importing the ATM waveform data
%                     w = [atm.shots(k).wf(i).w];
%
%   threshold         Amplitude threshold for centroid calculation. Threshold must be 0 < threshold <= 1. The ATM centroid tracker is using a 0.15 threshold
%                     of the maximum amplidude minus the noise floor. A threshold of 0.5 will outpute the full width at half maximum (FWHM) width of
%                     the pulse. 
%                     
% --------------------------------------------------------------------------------------------------------------------------------------------------------------
%
% OUTPUT PARAMETERS
%
%   t_ti              Time tags of the waveform in nano seconds that is equal or greater than the specified threhold. The first and last sample are 
%                     interpolated time tags between the regularly spaced t time tags at the specified threshold level.
%                     t_ti(end) - t_ti(1) is the width of the pulse at the specified threshold level. 
%
%   w_ti              Corresponding amplitude values for t_ti; 
%
%   t_c               Time tag of the centroid location in nano seconds.
%
% --------------------------------------------------------------------------------------------------------------------------------------------------------------
%
% Author:             Michael Studinger, NASA Goddard Space Flight Center, Greenbelt MD, USA.
% Version:            2.0 - December 18, 2020
% See also:           atm_wvfm_reader.m and atm_wvfm_info.m
% --------------------------------------------------------------------------------------------------------------------------------------------------------------

%% 
% turn everything into column vectors and convert to type double
% consider removing this check in future versions for faster execution

if iscolumn(w) == 1 % means is a column vector
    w = double(w);
elseif iscolumn(w) == 0
    w = double(w)';
end

if iscolumn(t) == 1 % means is a column vector
    t = double(t);
elseif iscolumn(t) == 0
    t = double(t)';
end

%% locate data above specified threshold

baseline = median(w(1:22)); % the noise floor is estimated as median of the first 21 smples.

[max_a, indx_max] = max(w);

% this fixes waveforms with negative values
if (min(w) < 0)
   w = w + abs(min(w)) + 1; 
end

w_t = threshold * (max_a - baseline) + baseline; 

% first start searching for first point which is above threshold going left
% from max amplitude
t_tmp = t(1:indx_max); % from first sample to maximum
w_tmp = w(1:indx_max);
% now reverse that
t_tmp_rev = wrev(t_tmp);
w_tmp_rev = wrev(w_tmp);
% now find first sample that's below threshold
indx_start_tmp = find((w_tmp_rev < w_t),1,'first'); 


% first start searching for first point which is above threshold going right from max amplitude
t_tmp = t(indx_max:end); % from first sample to maximum
w_tmp = w(indx_max:end);

% now find first sample that's below threshold 
indx_end_tmp = find((w_tmp < w_t),1,'first'); 

indx_start = indx_max - indx_start_tmp + 1;
indx_end = indx_max + indx_end_tmp - 1;

indx_above = indx_start:indx_end;

w_ti = w(indx_above);
t_i = t(indx_above);


%% do a linear interpolation around the first and last sample to estimate the entire curve above the specified threshold

t_s = ((w_t - w(indx_above(1))) .* (t(indx_above(2)) - t(indx_above(1))))./(w(indx_above(2)) - w(indx_above(1))) + t(indx_above(1));
t_e = ((w_t - w(indx_above(end-1))) .* (t(indx_above(end)) - t(indx_above(end-1))))./(w(indx_above(end)) - w(indx_above(end-1))) + t(indx_above(end-1));

w_ti(1) = w_t; w_ti(end) = w_t; % threshold per definition
t_i(1) = t_s;
t_i(end) = t_e;


%% use trapezoidal method to determine the approximate integral 

upper = trapz(t_i,(w_ti).*t_i);
lower = trapz(t_i,(w_ti));

t_c = upper/lower; % centroid

end


