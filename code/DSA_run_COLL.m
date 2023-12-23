function DSA = DSA_run_COLL(DSA)
%% DSA = DSA_run_COLL(DSA)
%
% This function initialises and runs the collocation method to determine
% a steady-state periodic (or constant) trajectory of the given nonlinear
% system in EQNS, for the set of numerical parameters defined in DATA.
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
    
    DSA = DSA_set(DSA);
    
    %% ---------------------------------------------------------------------
    % options:
    
    userJacobian    = strcmp(DSA.COLL.jacobian_type, 'analytical');
    algo_tol        = DSA.COLL.algo_tol;
    max_iter        = DSA.COLL.max_iter;

    userFunctionTolerance   = algo_tol;
    userStepTolerance       = algo_tol;
    userOptimalityTolerance = algo_tol;

    fsolve_opts = optimoptions('fsolve',...
        'FunctionTolerance',userFunctionTolerance, ...
        'StepTolerance',userStepTolerance, ...
        'OptimalityTolerance',userOptimalityTolerance, ...
        'MaxIterations',max_iter, ...
        'MaxFunctionEvaluations',1e9, ...
        'Algorithm','Levenberg-Marquardt', ...
        'Display','iter', ...
        'SpecifyObjectiveGradient',userJacobian);

    % ---------------------------------------------------------------------
    % (1) GET ANALYTICAL JACOBIAN
    if strcmp(DSA.COLL.jacobian_type, 'analytical')
        
        % Considered systems:
        % [dxdt] = [f(x,u,w)]
        % [z   ]   [g(x,u,w)]
        %  w     = z(t-T_d)
        %
        % J_RHS corresponds to the derivatives of [f(x,u,w); g(x,u,w)] with
        % respect to [x; w], with z the delayed variables and w the
        % to-delay variables.
        % J_LHS corresponds to the time derivation and delay matrices.
        
        % -- preparing J_RHS part, symbolically:
        if DSA.COLL.jacobian_reset
            DSA = DSA_get_NLSS(DSA); % get symbolic nonlinear system
            DSA_make_COLL_JAC(DSA); % get analytical jacobian
        end
        
        % -- preparing J_LHS part, numerically since it does not depend on
        % the analytical expressions of the system equations:
        
        % J_LHS_cell must contain N times matrix TTM_diff, followed by the 
        % D appropriate inverse-delay matrices. The code pays automatically
        % attention to the case of distinct delay values among delayed
        % variables.
        
        J_LHS_cell = {};
        J_LHS_cell(1:DSA.info.N) = {DSA.calc.TTM_diff};
        if DSA.data.delay.is_delayed
            for idx = 1:DSA.info.D
                this_delay = DSA.data.delay.(DSA.vars.lasi_todelay{idx});
                if this_delay == 0
                    J_LHS_cell(DSA.info.N + idx) = {eye(DSA.calc.nh_m)};
                else
                    % using "minus the delay" to obtain the inverse TTM_delay matrix:
                    J_LHS_cell(DSA.info.N + idx) = {DSA.calc.TTM_delay(-this_delay)};
                end
            end
        end
        DSA.calc.J_LHS = blkdiag(J_LHS_cell{:});
    end
    
    % ---------------------------------------------------------------------
    % (2) PREPARATION AND INITIALISATION OF HD DATA
    
    DSA = DSA_COLL_HD_init(DSA); % initialisation of HD data
    
    HD_COLL_1P = DSA.HD.HD_COLL_1P;
    M         = DSA.info.M;
    N         = DSA.info.N;
    D         = DSA.info.D;
    nh_m      = DSA.calc.nh_m;
    
    % ---------------------------------------------------------------------
    % (3) INITIALISATION OF UNKNOWN VECTOR
    % combining the init values of the unknown state and delay vars in a single vector
    
    % first initialising everything to zero:
    if DSA.data.delay.is_delayed
        varnames = [DSA.vars.lasi_states;
                    DSA.vars.lasi_delayed];
        unknown_vars_init = zeros((N+D)*nh_m, 1);
    else
        varnames = DSA.vars.lasi_states;
        unknown_vars_init = zeros(N*nh_m, 1);
    end
    
    % then filling in the initial guess:
    sz = size(varnames, 1);
    for vi = 1:sz
        vname = varnames{vi};
        if isfield(HD_COLL_1P, vname)
            tempval = real(DSA.calc.HD2TD*HD_COLL_1P.(vname));
        else
            tempval = 0.01; % any value should work, but preferably nonzero.
        end
        unknown_vars_init((vi-1)*nh_m + (1:nh_m)) = tempval;
    end
    
    % ---------------------------------------------------------------------
    % (4) INITIALISATION OF INPUT VECTOR
    % same as above but for the inputs
    % (This vector is called by the "get_u()" function.)
    
    varnames = DSA.vars.lasi_inputs;
    input_vars_vec = zeros(M*nh_m, 1);
    sz = size(varnames, 1);
    for vi = 1:sz
        vname = varnames{vi, 1};
        tempval = real(DSA.calc.HD2TD*HD_COLL_1P.(vname));
        input_vars_vec((vi-1)*nh_m + (1:nh_m)) = tempval;
    end
    DSA.COLL.input_vars_vec = input_vars_vec;
    
    % ---------------------------------------------------------------------
    % (5) COLLOCATION SYSTEM RESOLUTION
    % solving the nonlinear system
    
    DSA.calc.extract_all = 0;
    f_COLL = @(x) DSA_get_COLL_sstm(x, DSA);
    
    disp('--- COLL: start COLL solver (fsolve with Levenberg-Marquardt)')
    tic
    [sol_COLL_TD,FVAL,EXITFLAG,~,J_COLL] = fsolve(f_COLL, unknown_vars_init, fsolve_opts);
    toc
    
    if EXITFLAG > 0
        disp('--- COLL: solver is done')
    else
        warning('--- COLL: (possibly) NO SOLUTION FOUND. Check equations, data, or try another solver/algorithm.')
    end

    % ---------------------------------------------------------------------
    % (6) DISPLAY RESULTING ERROR

    figure
    hold on
    grid on
    stem(abs(FVAL))
    xlabel('unknowns')
    ylabel('residual error with respect to zero')
    title('Absolute residual COLL error')

    % ---------------------------------------------------------------------
    % (7) EXTRACTION
    % The COLL function is evaluated once again to obtain the spectrum of all
    % the variables in steady-state. Note that the initial conditions are
    % overwritten.
    
    disp('--- COLL: extracting results and waveforms')

    DSA.calc.extract_all = 1;
    HD_COLL_1P = DSA_get_COLL_sstm(sol_COLL_TD, DSA);
    
    % ---------------------------------------------------------------------
    % (8) TO TIME DOMAIN (resampling of the solution)
    
    TD_COLL_1P.t  = DSA.calc.DISP.t_1P;
    
    varnames = DSA.vars.lasi_all;
    for vi = 1:length(varnames)
        vname = varnames{vi};
        TD_COLL_1P.(vname)  = DSA_HD2TD(HD_COLL_1P.(vname), TD_COLL_1P.t,  HD_COLL_1P.f_b);
    end
    
    % ---------------------------------------------------------------------
    % (9) OUTPUTS
    % outputs and sorting structures
    DSA.HD.HD_COLL_1P  = HD_COLL_1P;
    DSA.TD.TD_COLL_1P  = TD_COLL_1P;
    DSA.COLL.J_COLL    = J_COLL;
    
    DSA.TD = orderfields(DSA.TD);
    DSA.HD = orderfields(DSA.HD);
    % -------------------------------------------------------------------------
    
    % if the system does not have pure delays, then the eigenvalues of the
    % jacobian matrix are the Floquet exponents
    if ~DSA.data.delay.is_delayed && DSA.COLL.jacobian_eigs
        EIGS_J_COLL = eig(J_COLL);
        
        disp('--- COLL: no non-zero delays: Jacobian eigenvalues are Floquet exponents')
        figure
        hold on
        grid on
        plot(real(EIGS_J_COLL), imag(EIGS_J_COLL)/(2*pi), '.', 'markersize', 12)
        if strcmp(get(groot, 'defaultTextInterpreter'), 'latex')
            xlabel('real part: $\sigma$ [rad/s]', 'fontsize', 14)
            ylabel('imag part: $\omega/(2\pi)$ [Hz]', 'fontsize', 14)
        else
            xlabel('real part: \sigma [rad/s]', 'fontsize', 14)
            ylabel('imag part: \omega/(2\pi) [Hz]', 'fontsize', 14)
        end
        yl = ylim;
        if abs(yl(1)-yl(2)) < 5e-4
           ylim(yl(1) + [-0.5 +0.5])
        end
        xl = xlim;
        if abs(xl(1)-xl(2)) < 5e-4
           xlim(xl(1) + [-0.5 +0.5])
        end
        title('eigenvalues of the jacobian')
    end
end

