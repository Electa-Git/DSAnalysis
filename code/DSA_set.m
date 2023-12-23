function DSA = DSA_set(DSA)
%% DSA = DSA_set(DSA)
% Setting up a dynamic system:
%
% - initialising data (DATA) and variables (VARS) from corresponding files
% - determination of system size
% - preparation of parameters for collocation method (COLL)
% - preparation of parameters for time-domain simulations (TDSIM)
%
% Author: Philippe De Rua
% Updated: 2023-06-06

%% Copyright © 2023 KU Leuven
% This file is part of DSAnalysis.
% DSAnalysis is a free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% DSAnalysis is distributed in the hope that it will be useful, but without 
% any warranty; without even the implied warranty of merchantability or 
% fitness for a particular purpose. See the GNU General Public License for 
% more details.
% You should have received a copy of the GNU General Public License along 
% with DSAnalysis. If not, see https://www.gnu.org/licenses/
    
% -------------------------------------------------------------------------
% initialising data and variables

DSA = DSA.fcts.DATA(DSA);
DSA = DSA.fcts.VARS(DSA);

% extracting fields:
vars  = DSA.vars;
info  = DSA.info;
calc  = DSA.calc;
fcts  = DSA.fcts;
COLL  = DSA.COLL;
TDSIM = DSA.TDSIM;

% It is assumed that the sets of states and delayed variables remains the 
% same before and after linearisation. Note that this may not be the case 
% in general, but this is a conservative approach: at least, no states nor
% delayed variables are lost in the process.

vars.smsi_states  = vars.lasi_states;
vars.smsi_delayed = vars.lasi_delayed;
vars.smsi_todelay = vars.lasi_todelay;

% -------------------------------------------------------------------------
% description of system size:
% 
% n - number of state variables
% m - number of input variables
% p - number of output variables
% d - number of delayed variables
%
% uppercase letters (N,M,P,D): lasi = large signal (nonlinear)  system
% lowercase letters (n,m,p,d): smsi = small signal (linearised) system

info.N = size(vars.lasi_states,  1);
info.M = size(vars.lasi_inputs,  1);
info.P = size(vars.lasi_outputs, 1);
info.D = size(vars.lasi_todelay, 1);

info.n = info.N;
info.m = size(vars.smsi_inputs,  1);
info.p = size(vars.smsi_outputs, 1);
info.d = info.D;

% safety checks:
if size(vars.lasi_todelay,1) ~= size(vars.lasi_delayed,1)
    error('There should be an equal number of todelay and delayed variables in VARS lists.')
end

% -------------------------------------------------------------------------
% automatic determination of unique delay values in system data:

% relevant fields are created:
%   delay.T:              a vector of all delay values
%   delay.unique_nonzero: a vector of all unique and strictly positive
%                         delay values. Duplicated delay values are removed.
%   delay.is_delayed:     a flag telling whether the system has any
%                         strictly positive delay.

if info.D > 0
    if isfield(DSA.data, 'delay')
        delay = DSA.data.delay;
    else
        error(['DSA.data.delay structure must be defined by the user in DATA' ...
            ' function. Every to-delay variable must be assigned with the numerical' ...
            ' value of its delay.'])
    end

    delay.T_d = zeros(1, info.D);

    % looping on the todelay variables:
    vnames = DSA.vars.lasi_todelay;
    for vi = 1:length(vnames)
        delay.T_d(vi) = delay.(vnames{vi});
    end
    
    % finding unique non-zero delay values:
    
    unique_T_d = unique(delay.T_d); % discarding duplicate delay values
    delay.unique_T_d = unique_T_d;
    
    % safety check:
    if any(delay.unique_T_d < 0)
        error('Delay values should be positive or zero.')
    end

    delay.unique_nonzero_T_d = unique_T_d(unique_T_d > 0); % discarding delays equal to zero
    delay.is_delayed = ~isempty(delay.unique_nonzero_T_d);
    delay.T_d_max = max(delay.unique_nonzero_T_d);
    delay.T_d_min = min(delay.unique_nonzero_T_d);
else
    % delay.unique_nonzero_T_d = []; not used? to be removed.
    delay.is_delayed = 0;
end

DSA.data.delay = delay;

% -------------------------------------------------------------------------
% size of harmonic vectors (COLL)

h_m = COLL.h_m;

if h_m < 0
    error('Parameter h_m cannot be negative. Please set a positive value (h_m >= 0). h_m = 0 will result in identifying a constant operating point. h_m > 0 results in identifying a possibly periodic trajectory.')
end

nh_m      = 2*h_m+1; % total number of components
calc.h_m  = h_m;
calc.nh_m = nh_m;

% -------------------------------------------------------------------------
% pre-calculations related to time and frequency data:

f_b     = DSA.data.f_1; % base frequency is by default equal to system fundamental frequency [Hz]
omega_b = 2*pi*f_b;     % base angular frequency [rad/s]
T_b     = 1/f_b;        % fundamental period [s]

calc.h       = -h_m:h_m; % vector of harmonic indices
calc.f       = f_b * calc.h; % corresponding vector of frequency values
calc.f_b     = f_b;
calc.omega_b = omega_b;
calc.T_b     = T_b;

% --- collocation method calculations (COLL)

COLL.T_s = T_b / nh_m; % sampling time [s]
COLL.f_s = nh_m * f_b; % sampling frequency [Hz]

% --- resampled collocation waveforms (COLL)

L_DISP    = COLL.L_DISP;  % number of sampling points per fundamental period
DISP.T_s  = T_b / L_DISP; % sampling period [s]
DISP.f_s  = L_DISP * f_b; % sampling frequency [Hz]

% --- time domain simulations/numerical integrations (TDSIM)

tspan     = TDSIM.tspan;
XP        = TDSIM.XP;
tduration = tspan(end) - tspan(1);
NP = tduration/T_b;
if XP > NP
    XP = floor(NP);
    if XP == 0
        XP = 1;
    end
end

TDSIM.XP  = XP;
TDSIM.NP  = NP;
TDSIM.f_s = 1/TDSIM.T_s; % sampling frequency [Hz]

% --- definition of time vectors [s]
% 
% 1P     - one fundamental period, [0,    T_b[, last point excluded
% XP     - XP fundamental periods, [0, XP*T_b[, last point excluded
% NP     - NP fundamental periods, [0, NP*T_b], last point included

step       = COLL.T_s;
stop       = T_b - COLL.T_s;
COLL.t_1P  = (0:step:stop)';  

step       = T_b/L_DISP;
stop       = T_b - T_b/L_DISP;
DISP.t_1P  = (0:step:stop)';

step       = TDSIM.T_s;
stop       = T_b - TDSIM.T_s;
TDSIM.t_1P = (0:step:stop)';

step       = TDSIM.T_s;
stop       = XP*T_b - TDSIM.T_s;
TDSIM.t_XP = (0:step:stop)';

step       = TDSIM.T_s;
stop       = NP*T_b;
TDSIM.t_NP = tspan(1) + (0:step:stop)';

% -------------------------------------------------------------------------
% calculation of time- and frequency-transfer matrices for the COLL method

% shifted, i.e. for harmonic vectors with spectrum [0 ... hm -hm ... -1]
TD2HD_shifted = dftmtx(nh_m)/nh_m;
HD2TD_shifted = conj(dftmtx(nh_m));
% centered, i.e. for harmonic vectors with spectrum [-hm ... -1 0 ... hm]
TD2HD = TD2HD_shifted([h_m+2:nh_m, 1:h_m+1],:); 
HD2TD = HD2TD_shifted(:, [h_m+2:nh_m, 1:h_m+1]);

% Frequency transfer matrix for derivatives and for delays (i.e. just a phase shift in Harmonic Domain)
FTM_diff  = 1j*omega_b * diag(-h_m:h_m);
FTM_delay = @(T_d) diag(exp(-1j*omega_b*(-h_m:h_m)*T_d));

% Time transfer matrices for derivatives and for delays
% the resulting matrices *are* real, so just use real() to remove the 
% remaining complex part caused by machine precision.
TTM_diff  =        real(HD2TD * FTM_diff * TD2HD);       % the result *is* real, so just remove the complex part caused by machine precision.
TTM_delay = @(T_d) real(HD2TD * FTM_delay(T_d) * TD2HD); % the result *is* real, so just remove the complex part caused by machine precision.

calc.TTM_diff  = TTM_diff;
calc.TTM_delay = TTM_delay;
calc.TD2HD     = TD2HD;
calc.HD2TD     = HD2TD;

% -------------------------------------------------------------------------
% referencing of variable indices:
vars.lasi_all = [
     vars.lasi_states;
     vars.lasi_inputs;
     vars.lasi_outputs;
     vars.lasi_todelay;
     vars.lasi_delayed;
     ];

varsets = {
 'lasi_states';
 'lasi_inputs';
 'lasi_outputs';
 'lasi_todelay';
 'lasi_delayed';
 'smsi_inputs';
 'smsi_outputs';
};

for iset = 1:size(varsets,1)
    varset = varsets{iset};
    vnames = vars.(varset);
    for vi = 1:size(vnames, 1)
        vname = vnames{vi, 1};
        lidx.(varset).(vname) = vi;
    end
end

vars.lidx = lidx;

% -------------------------------------------------------------------------
% mapping functions:

% mapping from harmonic to matlab index:
fcts.h2i_h_m = @(h) h + h_m + 1;

% finding unknown names based on unknown index:
%vars.lasi_unknowns = [vars.lasi_states; vars.lasi_delayed];
lasi_unknowns = [vars.lasi_states; vars.lasi_delayed];
fcts.i2unknown = @(i) lasi_unknowns{floor((i-1)/calc.nh_m)+1};

% -------------------------------------------------------------------------
% gathering output:

DSA.calc      = calc;
DSA.info      = info;
DSA.vars      = vars;
DSA.fcts      = fcts;
DSA.COLL      = COLL;
DSA.TDSIM     = TDSIM;
DSA.calc.DISP = DISP;

fprintf('\n--- DSA_set: DSA updated\n')
end