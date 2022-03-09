function f = atm_gaussian_two_peak(pars,t)
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% SUMMARY:          The atm_gaussian_two_peak function can be used as a function or a function handle for input into a nonlinear regression.
%                   
% SYNTAX:           use as a function handle: @atm_gaussian_two_peak
%                   use as a function:        atm_gaussian_two_peak(pars,t)
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% INPUT PARAMETERS
% 
%  pars             Array with 7 input parameters for a Gaussian function with two peaks as follows: 
%                   e_g1 = pars(1); % signal baseline 1
%                   a_g1 = pars(2); % amplitude 1
%                   t_g1 = pars(3); % peak location 1
%                   s_g1 = pars(4); % 1-sigma peak width = standard deviation
% 
%                   a_g2 = pars(5); % amplitude 2
%                   t_g2 = pars(6); % peak location 2
%                   s_g2 = pars(7); % 1-sigma peak width = standard deviation
%                   
%  t                Array of ATM waveform time tags with laser trigger time in nano seconds (e.g., atm_wvfm.shots(shot_nr).wf(rcv_gate).t)  
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% OUTPUT PARAMETERS 
%                   function f when used as function handle or
%                   result of nonlinear regression when used as a function
% 
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Author:           Michael Studinger, NASA Goddard Space Flight Center, Greenbelt MD, USA.
% Version:          1.02 - August 19, 2015
% See also:         
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------

%% check input parameters

% check number of input arguments
if (nargin ~= 2)    
        error('atm_gaussian_two_peak:nargin', ['\n\tERROR: Number of input arguments must be 2:\n', ...
        '\tSYNTAX: atm_gaussian_two_peak(pars,t);\n'])
end

% check pars array
if (any(~isfinite(pars)) || ~isvector(pars) || numel(pars) ~= 7)
    error('atm_gaussian_two_peak:pars', '\n\tERROR: pars input values must be finite and a 7 element vector.\n');
end

if (isvector(t) == 1)     % means is a not column vector (dimension 1 x n)
    t = double(t);
elseif (isvector(t) == 0) % is a column vector (dimension n x 1) that needs to be transformed into a vector 
    t = double(t)';
end

%% set Gaussian function with two peaks and one baseline

e_g1 = pars(1); % signal baseline 1
a_g1 = pars(2); % amplitude 1
t_g1 = pars(3); % peak location 1
s_g1 = pars(4); % 1-sigma peak width = standard deviation

a_g2 = pars(5); % amplitude 2
t_g2 = pars(6); % peak location 2
s_g2 = pars(7); % 1-sigma peak width = standard deviation

f =  e_g1 + a_g1.* exp(-0.5.*(t - t_g1).^2/(s_g1^2)) + a_g2.* exp(-0.5.*(t - t_g2).^2/(s_g2^2));

end

