function [output, J] = DSA_get_COLL_sstm(x_unknowns, DSA)
%% [output, J] = DSA_get_COLL_sstm(x_unknowns, DSA)
% 
% Interface function to obtain the collocation formulation from the EQNS
% file.
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
    
    %#ok<*NASGU>

    % ---------------------------------------------------------------------
    % extraction:
    
    delay      = DSA.data.delay;
    is_delayed = DSA.data.delay.is_delayed;
    nh_m       = DSA.calc.nh_m;
    TTM_diff   = DSA.calc.TTM_diff;
    TTM_delay  = DSA.calc.TTM_delay;
    
    % extraction of parameters :
    param = DSA.data.param;

    % extraction of state and delayed variables:
    get_x = @(v) x_unknowns((v-1)*nh_m + (1:nh_m)); % extracts as vector
    
    % extraction of inputs:
    get_u = @(v) DSA.COLL.input_vars_vec((v-1)*nh_m + (1:nh_m));
    
    % -------------------------------------------------------------------------
    % calculations:
    
    get_values  = 1; % --> need numerical values of parameters and variables
    run_algebra = 1; % --> running algebraic relationships leading to f(x,u)
    run_dxdt    = ~DSA.calc.extract_all; % --> getting the expression of dxdt = f(x,u), except for the last run
    
    % function get_u fetches input variables from their index:
    for idx = 1:DSA.info.M
        inputs.(DSA.vars.lasi_inputs{idx}) = get_u(idx); %#ok<STRNU> 
    end
    % function get_x fetches unknown variables from their index:
    for idx = 1:DSA.info.N
        states.(DSA.vars.lasi_states{idx}) = get_x(idx); %#ok<STRNU> 
    end
    for idx = 1:DSA.info.D
        if is_delayed
            delayed.(DSA.vars.lasi_delayed{idx}) = get_x(DSA.info.N+idx);
        else
            delayed.(DSA.vars.lasi_delayed{idx}) = NaN;
        end
    end
    
    run(DSA.fcts.EQNS)
    
    % -------------------------------------------------------------------------

    if DSA.calc.extract_all
        % (this is run after the COLL solver is done)

        HD_COLL_1P = DSA.HD.HD_COLL_1P;
        TD2HD = DSA.calc.TD2HD;
        
        % getting the list of all variables to export:
        if is_delayed
            varnames = [DSA.vars.lasi_states;
                        DSA.vars.lasi_inputs;
                        DSA.vars.lasi_outputs;
                        DSA.vars.lasi_todelay;
                        DSA.vars.lasi_delayed];
        else
            varnames = [DSA.vars.lasi_states;
                        DSA.vars.lasi_inputs;
                        DSA.vars.lasi_outputs];
        end
        for vi = 1:length(varnames)
            vname = varnames{vi, 1};
            % checking whether the variables exist:
            if ~exist(vname, 'var')
                warning(['Variable ' vname ' is not defined in EQNS. This may cause an error. Either remove the variable from the list in VARS, or provide its analytical expression in EQNS.'])
                continue
            end
            % Bringing the solution from the (time-domain) COLL solution to
            % the HD:
            TD_vec = eval(vname);
            HD_COLL_1P.(vname) = TD2HD*TD_vec;
        end
        % If the system allows for pure delays but all delays are set to
        % zero, then the delayed variables are not defined. The solution to 
        % this problem is to set them equal to their todelay counterparts:
        if ~is_delayed && DSA.info.D > 0
            for vi = 1:DSA.info.D
                vname_delayed = DSA.vars.lasi_delayed{vi};
                vname_todelay = DSA.vars.lasi_todelay{vi};
                TD_vec_todly = eval(vname_todelay);
                HD_COLL_1P.(vname_todelay) = TD2HD*TD_vec_todly;
                HD_COLL_1P.(vname_delayed) = HD_COLL_1P.(vname_todelay);
            end
        end
        
        output = HD_COLL_1P;
        J = NaN; % (The last-iteration jacobian is given by fsolve)
        
    else
        % This is run while COLL solver is working, i.e. at every iteration
        % step.
        
        % Initialisation of left-hand side (LHS) and right-hand side (RHS)
        % of dxdt = f(x,u):
        
        LHS_dxdt = zeros(size(x_unknowns));
        RHS_dxdt = LHS_dxdt; % (just to get the size right)
        RHS_dxdt(1:(DSA.info.N*nh_m)) = dxdt;

        % calculation of the time-derivative of the state variables:
        for idx = 1:DSA.info.N
            LHS_dxdt((1:nh_m) + (idx-1)*nh_m) = TTM_diff * get_x(idx);
        end

        if DSA.data.delay.is_delayed
            % extending with the delayed and non-delayed versions of delay variables:
            for idx = 1:DSA.info.D
                this_delay = DSA.data.delay.(DSA.vars.lasi_todelay{idx});
                if this_delay == 0
                    LHS_dxdt((1:nh_m) + (DSA.info.N+idx-1)*nh_m) = get_x(DSA.info.N+idx);
                else
                    % using "minus the delay" to obtain the inverse TTM_delay matrix:
                    LHS_dxdt((1:nh_m) + (DSA.info.N+idx-1)*nh_m) = TTM_delay(-this_delay) * get_x(DSA.info.N+idx);
                end
                RHS_dxdt((1:nh_m) + (DSA.info.N+idx-1)*nh_m) = eval(DSA.vars.lasi_todelay{idx});
            end
        end
        
        delta = RHS_dxdt - LHS_dxdt; % the solution error being minimised by fsolve
        output = delta;
        
        if nargout == 2
            % the solver is calling for the jacobian.
            % Jacobian calculation: there are two terms: J = J_RHS - J_LHS
            %
            % J_LHS comes from the derivatives of periodic variables N*X
            % J_LHS is precalculated numerically and stored in DSA.calc
            %
            % J_RHS comes from F(X,U)
            % J_RHS evaluated numerically with the following runs:
            
            null_symb = zeros(nh_m,1);
            if DSA.data.delay.is_delayed
                run(DSA.fcts.COLL_JAC_DLY)
            else
                run(DSA.fcts.COLL_JAC)
            end
            
            J = real(J_RHS - DSA.calc.J_LHS);
        end
    end
    
    % -------------------------------------------------------------------------
end
   