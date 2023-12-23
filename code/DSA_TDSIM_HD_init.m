function DSA = DSA_TDSIM_HD_init(DSA)
%% DSA = DSA_TDSIM_HD_init(DSA)
%
% This function selects the harmonic spectrum of states and delayed 
% variables when initialising the numerical integration (TDSIM) from a a
% periodic trajectory, i.e. when initialising from either 'source' or
% 'user', which are harmonic domain definitions of periodic trajectories
% (their spectra).
%
% Author: Philippe De Rua
% Updated: 2023-11-28

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

f_b  = DSA.calc.f_b;
vars = DSA.vars;
h_m  = DSA.calc.h_m;
nh_m = DSA.calc.nh_m;

% ---------------------------------------------------------------------
% (0.1) CREATION OF EMPTY STRUCTURE

varnames = [vars.lasi_states;
            vars.lasi_delayed];
for vi = 1:size(varnames, 1)
    vname = varnames{vi, 1};
    empty_HD.(vname) = zeros(nh_m, 1); % creation and init to zero
end

% ---------------------------------------------------------------------
% (0.2) INITIALISING UNKNOWNS

if strcmp(DSA.TDSIM.init_from, 'source')
    % GETTING SOURCE DATA (DSA.TDSIM.HD_source)

    if isfield(DSA, 'HD') && isfield(DSA.HD, DSA.TDSIM.HD_source)
        src_HD = DSA.HD.(DSA.TDSIM.HD_source);
        
        % checking that the base frequency used for all calculation is
        % compatible with the user's data:
        if src_HD.f_b ~= f_b
            error('--- DSA_TDSIM_HD_init: data.f_b is not compatible with HD_source data')
        end
    
    elseif strcmp(DSA.TDSIM.init_from, 'source')
        error(['--- DSA_TDSIM_HD_init: HD_source ' DSA.TDSIM.HD_source ' is not available in DSA.HD structure'])
    end

    disp('--- DSA_TDSIM_HD_init: initialising unknowns from HD source:')
    old_HD = src_HD;

elseif strcmp(DSA.TDSIM.init_from, 'user')
    % GETTING USER DATA (HD)

    HDi = empty_HD; % (giving it empty)
    usr_HD = DSA.fcts.HD(DSA, HDi); % (gitting it back, filled)
    
    % checking that the base frequency used for all calculation is
    % compatible with the user's data:
    if usr_HD.f_b ~= f_b
        error('--- DSA_TDSIM_HD_init: data.f_b is not compatible with user''s HD data')
    end

    disp('--- DSA_TDSIM_HD_init: initialising unknowns from user''s HD function:')
    old_HD = usr_HD;

else
    error(['--- DSA_TDSIM_HD_init: cannot initialise unknowns from ' DSA.TDSIM.init_from])
end

varnames = [vars.lasi_states;
            vars.lasi_delayed];
new_HD = empty_HD;
new_HD = DSA_HD_init_from_source(old_HD, new_HD, varnames, h_m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in case h_m has changed:
new_HD.f    = f_b*(-h_m:h_m);
new_HD.h_m  = h_m;
new_HD.nh_m = nh_m;
new_HD.f_b  = f_b;

% returning the new_HD
DSA.TDSIM.HD_TDSIM_1P_init = new_HD;

disp('--- DSA_TDSIM_HD_init: finished HD data initialisation')
end