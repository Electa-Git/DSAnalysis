function DSA = DSA_COLL_HD_init(DSA)
%% DSA = DSA_COLL_HD_init(DSA)
%
% This function initialises the harmonic-domain data, i.e. the harmonic
% spectrum of the system variables. The data is initialised from the
% user-defined data in HD file, from the latest solution of the COLL method
% or from any other HD structure provided by the user.
% 
% The inputs have a predefined harmonic content. The harmonic spectrum of
% states and delayed variables is a priori unknown and will be determined
% with the collocation (COLL) method. Consequently, the spectrum of these
% unknown variables is initialised as an initial guess for the solver.
% 
% When initialising the HD, attention is paid to the maximum harmonic rank
% h_m requrested from the collocation method. Initialising from a source
% with smaller h_m is safe. Initialising from a source with larger h_m
% means that information from the source is lost.
%
% Variables defined in the lists (VARS) but without numerical data in HD
% file or in the user-provided HD source are set to zero by default. Pay
% attention to the fact that the variables names are case sensitive.
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

f_b  = DSA.calc.f_b;
vars = DSA.vars;
h_m  = DSA.calc.h_m;
nh_m = DSA.calc.nh_m;

% ---------------------------------------------------------------------
% (0.0) CREATION OF EMPTY STRUCTURE

varnames = [vars.lasi_inputs;
            vars.lasi_states;
            vars.lasi_delayed];
for vi = 1:size(varnames, 1)
    vname = varnames{vi, 1};
    empty_HD.(vname) = zeros(nh_m, 1); % creation and init to zero
end

% ---------------------------------------------------------------------
% (0.1) SAVING LATEST SOLUTION

if isfield(DSA.HD, 'HD_COLL_1P')
    lastsol_HD = DSA.HD.HD_COLL_1P;
elseif strcmp(DSA.COLL.init_guess_from, 'solution') || strcmp(DSA.COLL.def_inputs_from, 'solution')
    warning('DSA_COLL_HD_init: HD_COLL_1P is not available in DSA.HD, so an initialisation from the last solution is not possible. Try setting ''init_guess_from'' to ''user'' instead. For this particular run, the code will initialise the harmonic spectra from the user-defined HD file.')
    DSA.COLL.init_guess_from = 'user';
end

% ---------------------------------------------------------------------
% (0.2) GETTING USER'S DATA

HDi = empty_HD; % (giving it empty)
usr_HD = DSA.fcts.HD(DSA, HDi); % (getting it back, filled)

% ---------------------------------------------------------------------
% (0.3) GETTING SOURCE DATA (DSA.COLL.HD_source)

if isfield(DSA.HD, DSA.COLL.HD_source)
    src_HD = DSA.HD.(DSA.COLL.HD_source);
    
elseif strcmp(DSA.COLL.init_guess_from, 'source') || strcmp(DSA.COLL.def_inputs_from, 'source')
    error(['--- DSA_COLL_HD_init: HD_source ' DSA.COLL.HD_source ' is not available in DSA.HD structure'])
end

% ---------------------------------------------------------------------
% (0.4) INITIALISING UNKNOWNS 
% accounting for possibly different h_m values:

if     strcmp(DSA.COLL.init_guess_from, 'solution')
    
    disp('--- DSA_COLL_HD_init: initialising unknowns from last solution:')
    old_HD = lastsol_HD;

elseif strcmp(DSA.COLL.init_guess_from, 'source')
    
    % checking that the base frequency used for all calculation is
    % compatible with the user's data:
    if src_HD.f_b ~= f_b
        error('--- DSA_COLL_HD_init: base frequency DSA.calc.f_b is not compatible with HD_source data')
    end
    
    disp('--- DSA_COLL_HD_init: initialising unknowns from HD source:')
    old_HD = src_HD;

elseif strcmp(DSA.COLL.init_guess_from, 'user')

    % checking that the base frequency used for all calculation is
    % compatible with the user's data:
    if usr_HD.f_b ~= f_b
        error('--- DSA_COLL_HD_init: base frequency DSA.calc.f_b is not compatible with user''s HD data')
    end

    disp('--- DSA_COLL_HD_init: initialising unknowns from user-defined HD function:')
    old_HD = usr_HD;

else
    error(['--- DSA_COLL_HD_init: cannot initialise unknowns from ' DSA.COLL.init_guess_from])
end

varnames =  [vars.lasi_states;
             vars.lasi_delayed];
new_HD = empty_HD;
new_HD = DSA_HD_init_from_source(old_HD, new_HD, varnames, h_m);

% ---------------------------------------------------------------------
% (0.5) INITIALISING INPUTS 
% accounting for possibly different h_m values:
    
if     strcmp(DSA.COLL.def_inputs_from, 'solution')
    
    disp('--- DSA_COLL_HD_init: defining inputs from last solution:')
    old_HD = lastsol_HD;

elseif strcmp(DSA.COLL.def_inputs_from, 'source')
    
    disp('--- DSA_COLL_HD_init: defining inputs from HD source:')
    old_HD = src_HD;

elseif strcmp(DSA.COLL.def_inputs_from, 'user')
    
    disp('--- DSA_COLL_HD_init: defining inputs from user-defined HD function:')
    old_HD = usr_HD;

else
    error(['--- DSA_COLL_HD_init: cannot initialise inputs from ' DSA.COLL.def_inputs_from])
end

varnames = vars.lasi_inputs;
new_HD = DSA_HD_init_from_source(old_HD, new_HD, varnames, h_m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in case h_m has changed:
new_HD.f    = f_b*(-h_m:h_m);
new_HD.h_m  = h_m;
new_HD.nh_m = nh_m;
new_HD.f_b  = f_b;

% returning the new_HD as HD_COLL_1P
DSA.HD.HD_COLL_1P = new_HD;

disp('--- DSA_COLL_HD_init: finished HD data initialisation')
end