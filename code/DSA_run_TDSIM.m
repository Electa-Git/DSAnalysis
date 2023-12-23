function DSA = DSA_run_TDSIM(DSA)
%% DSA = DSA_run_TDSIM(DSA)
%
% This function performs a time-domain simulation of the nonlinear dynamic
% system via numerical integration of the differential equations, starting
% from initial conditions (ODEs) or an initial function segment (DDEs).
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
    
    DSA = DSA_set(DSA);

    %% ------------------------------------------------------------------------
    % extraction of necessary data
    % -------------------------------------------------------------------------
    
    f_s        = DSA.TDSIM.f_s;
    t_1P       = DSA.TDSIM.t_1P;
    t_XP       = DSA.TDSIM.t_XP;
    t_NP       = DSA.TDSIM.t_NP;
    XP         = DSA.TDSIM.XP;
    NP         = DSA.TDSIM.NP;
    h_m        = DSA.calc.h_m;
    delay      = DSA.data.delay;
    is_delayed = DSA.data.delay.is_delayed;
    
    DSA.TDSIM.t_start = 0; % for now, but the code can be improved to start
                           % other time values than 0.

    DSA.calc.extract_all = false; % will only be true once the simulation is complete
    sim_opts         = odeset('RelTol',DSA.TDSIM.rel_tol,'AbsTol',DSA.TDSIM.abs_tol,'OutputFcn',@DSA_odeprog);
    sim_opts_no_prog = odeset('RelTol',DSA.TDSIM.rel_tol,'AbsTol',DSA.TDSIM.abs_tol);

    %% ------------------------------------------------------------------------
    % automatic definition of functions
    % -------------------------------------------------------------------------

    if ~is_delayed
        fun_dxdt = @(t, x) DSA_dxdt_lasi(t, x, DSA);
    else
        fun_dxdt_delay = @(t, x, X_delay) DSA_dxdt_lasi_delay(t, x, X_delay, DSA);
    end

    %% ------------------------------------------------------------------------
    % preparation of the state variables initial conditions
    % -------------------------------------------------------------------------
    
    % The large-signal simulations start from a limit cycle or from zero :

    if strcmp(DSA.TDSIM.init_from, 'source') || strcmp(DSA.TDSIM.init_from, 'user')
        DSA = DSA_TDSIM_HD_init(DSA);
        HD_struct = DSA.TDSIM.HD_TDSIM_1P_init;
        
        if ~is_delayed
            x_init =      DSA_xOft_hist_lasi(0, DSA.vars.lasi_states, HD_struct);
        else
            x_hist = @(t) DSA_xOft_hist_lasi(t, DSA.vars.lasi_states, HD_struct);
        end

    elseif strcmp(DSA.TDSIM.init_from, 'solution')
        if ~is_delayed
            if isfield(DSA.TDSIM, 'x_solution')
                x_init = DSA.TDSIM.x_solution;
            else
                error('--- TDSIM: initialisation from "solution" is not possible: no solution available yet.')
            end
        else
            if isfield(DSA.TDSIM, 'x_solution')
                x_hist = DSA.TDSIM.x_solution;
            else
                error('--- TDSIM: initialisation from "solution" is not possible: no solution available yet.')
            end
        end

    elseif strcmp(DSA.TDSIM.init_from, 'zero')
        if ~is_delayed
            x_init = zeros(DSA.info.N, 1);
        else
            x_hist = zeros(DSA.info.N, 1);
        end

    else
        error(['--- TDSIM: initialisation from "' DSA.TDSIM.init_from '" is not available'])
    end

    disp('--- TDSIM: preparation done')

    %% ------------------------------------------------------------------------
    % solving the systems
    % -------------------------------------------------------------------------

    if ~is_delayed
        disp(['--- TDSIM: solver starts with ' DSA.TDSIM.solver])
        [~, x_lasi] = feval(DSA.TDSIM.solver, fun_dxdt, t_NP, x_init, sim_opts);

        DSA.TDSIM.x_solution = x_lasi(end,:).';

    else % solving with dde23
        disp('--- TDSIM: solver starts with dde23')
        T_d_min = delay.T_d_min;
        
        % the total time span is split into subintervals with length
        % equal to the smallest delay. The total number of subintervals is:
        nb_subint = ceil((t_NP(end)-t_NP(1))/T_d_min);
    
        % solve via dde23:
        % to avoid very long simulation times, it is better to split 
        % the integration into subintervals. dde23_sol accumulates
        % the solution info over all subintervals.
        
        % looping on the subintervals
        disp('--- TDSIM: integrating subinterval:')
        for subint_idx = 1:nb_subint
            % displaying progress:
            if (subint_idx == 1) || (subint_idx == nb_subint)
                fprintf('%5d on %5d\n', subint_idx, nb_subint)
            elseif mod(subint_idx, 10) == 0
                fprintf('%5d on %5d\n', subint_idx, nb_subint)
            end
            
            % determining time values of this subinterval:
            if subint_idx < nb_subint
                t_ddly_start = t_NP(1) + (subint_idx-1)*T_d_min;
                t_ddly_stop  = t_NP(1) + subint_idx*T_d_min;
            else
                t_ddly_start = t_NP(1) + (subint_idx-1)*T_d_min;
                t_ddly_stop  = t_NP(end);
            end
            t_ddly = unique([t_ddly_start; t_NP((t_ddly_start <= t_NP) & (t_NP <= t_ddly_stop)); t_ddly_stop]);
            
            % solving on the subinterval:
            if subint_idx == 1
                x_prev_sol = x_hist;
            else % subint_idx > 1
                x_prev_sol = dde23_sol;
            end
            dde23_sol = dde23(fun_dxdt_delay, delay.unique_nonzero_T_d, x_prev_sol, t_ddly, sim_opts_no_prog);
        end
        
        DSA.TDSIM.x_solution = dde23_sol;

        disp('--- TDSIM: extracting state variables')
        
        % extracting the non-delayed state variables:
        
        warning('off', 'MATLAB:deval:NonuniqueSolution')
        x_lasi = deval(dde23_sol, t_NP);
        warning('on', 'MATLAB:deval:NonuniqueSolution')

        % extracting the delayed state variables:

        nb_unique_dlys = length(delay.unique_nonzero_T_d);
        X_dly_lasi = zeros(DSA.info.N, nb_unique_dlys, length(t_NP)); 
        for dly_idx = 1:nb_unique_dlys % for each unique value of time delay

            % 1- determining the delay value:
            T_di = delay.unique_nonzero_T_d(dly_idx);

            % 2- complete time vector over which delayed variables have values:
            t_ddly_Tdi = t_NP - T_di;

            % 3- splitting time vector between negative and positive time values:
            neg_t_ode_Tdi = t_ddly_Tdi(t_ddly_Tdi <  t_NP(1));
            pos_t_ode_Tdi = t_ddly_Tdi(t_ddly_Tdi >= t_NP(1));

            % 4- processing negative time values with the history function:
            if isempty(neg_t_ode_Tdi)
                neg_part = [];
            else
                if isa(x_hist,'function_handle')
                    neg_part = x_hist(neg_t_ode_Tdi);
                elseif isa(x_hist,'struct')
                    warning('off', 'MATLAB:deval:NonuniqueSolution')
                    neg_part = deval(x_hist, neg_t_ode_Tdi);
                    warning('on', 'MATLAB:deval:NonuniqueSolution')
                else
                    % in this case, the vector must be repeated manually
                    neg_part = repmat(x_hist, 1, length(neg_t_ode_Tdi));
                end
            end

            % 5- processing positive time values with the deval function:
            warning('off', 'MATLAB:deval:NonuniqueSolution')
            pos_part = deval(dde23_sol, pos_t_ode_Tdi);
            warning('on', 'MATLAB:deval:NonuniqueSolution')

            % 6- combining positive and negative times:
            X_dly_lasi(:, dly_idx, :) = [neg_part pos_part];
        end
    end
    
    disp('--- TDSIM: solver is done')

    %% ------------------------------------------------------------------------
    % extracting input and output signals
    % -------------------------------------------------------------------------

    disp('--- TDSIM: extracting input and output variables')

    % Extracting of input and output variables:
    DSA.calc.extract_all = true;
    if ~is_delayed
        dxdt_output = DSA_dxdt_lasi(t_NP, x_lasi, DSA);
    else
        dxdt_output = DSA_dxdt_lasi_delay(t_NP, x_lasi, X_dly_lasi, DSA);
    end
    TD_TDSIM_NP = dxdt_output.TD_DE;
    TD_TDSIM_NP.t = t_NP;
    
    DSA.TD.TD_TDSIM_NP = TD_TDSIM_NP;

    % ---------------------------------------------------------------------
    % Extracting the last fundamental period, if t_extr is at least
    % as long as the fundamental period. Then, extracting the harmonic
    % content of the last fundamental period.
    if NP >= 1
        mask_1P = (length(t_NP)-length(t_1P)+1):(length(t_NP));
        if all(mask_1P > 0)
            disp('--- TDSIM: Fourier analysis of the last fundamental period')
            varnames = DSA.vars.lasi_all;
            for vi = 1:size(varnames)
                vname = varnames{vi};
                TD_TDSIM_1P.(vname) = TD_TDSIM_NP.(vname)(mask_1P);
            end
            TD_TDSIM_1P.t = t_1P;
        
            [HD_TDSIM_1P, FHD_TDSIM_1P] = DSA_TD2HD_4struct(TD_TDSIM_1P, varnames, f_s, h_m);
            DSA.TD.TD_TDSIM_1P  = TD_TDSIM_1P;
            DSA.HD.HD_TDSIM_1P  = HD_TDSIM_1P;
            DSA.HD.FHD_TDSIM_1P = FHD_TDSIM_1P;
            
            if DSA.HD.FHD_TDSIM_1P.f_b ~= DSA.calc.f_b
                warning('The sampling time should divide the fundamental period into an integer number of steps for the frequency increment df to be an integer divider of the fundamental frequency f_1.')
            end
        else
            warning('--- TDSIM: time range is not long enough to extract the last fundamental period and related harmonic content.')
        end
    else
        warning('--- TDSIM: time range is not long enough to extract the last fundamental period and related harmonic content.')
    end

    % ---------------------------------------------------------------------
    % Extracting the last XP fundamental periods and corresponding harmonic
    % content:
    if NP >= XP
        mask_XP = (length(t_NP)-length(t_XP)+1):(length(t_NP));
        if all(mask_XP > 0)
            disp('--- TDSIM: Fourier analysis of the last XP fundamental periods')
            varnames = DSA.vars.lasi_all;
            for vi = 1:size(varnames)
                vname = varnames{vi};
                TD_TDSIM_XP.(vname) = TD_TDSIM_NP.(vname)(mask_XP);
            end
            TD_TDSIM_XP.t = t_XP;
            
            [HD_TDSIM_XP, FHD_TDSIM_XP] = DSA_TD2HD_4struct(TD_TDSIM_XP, varnames, f_s, h_m);
            DSA.TD.TD_TDSIM_XP  = TD_TDSIM_XP;
            DSA.HD.HD_TDSIM_XP  = HD_TDSIM_XP;
            DSA.HD.FHD_TDSIM_XP = FHD_TDSIM_XP;
        else
            warning('--- TDSIM: time range is not long enough to extract the last XP fundamental period(s) and related harmonic content.')
        end
    else
        warning('--- TDSIM: time range is not long enough to extract the last XP fundamental period(s) and related harmonic content.')
    end

    % ---------------------------------------------------------------------
    % sorting structures
    DSA.TD = orderfields(DSA.TD);
    DSA.HD = orderfields(DSA.HD);

    disp('--- TDSIM: signals extraction is done')
end