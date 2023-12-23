function DSA = DSA_get_HSS(DSA)
%% DSA = DSA_get_HSS(DSA)
%
% Function DSA_get_HSS fills the DSA structure with field LTI (see DSA.LTI
% afterwards. This field contains the LTI system representation. Generally,
% this system corresponds to a frequency-lifted system, since it is the
% main purpose of this toolbox to deliver such systems. However, it is also
% possible to work with non-lifted LTI systems, by simply setting parameter
% h_t (the truncation rank) to zero. Please refer to the PhD thesis of the
% author to understand how to choose the truncation rank, as well as its
% relation to the periodic rank and the forced periodic rank.
% 
% This function is essentially a dispatch function, which calls either 
% DSA_LTP2HSS(DSA) or DSA_LTP2DHSS(DSA) according to whether the system has
% pure delays or not.
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

%#ok<*UNRCH>

    %% Transforming the LTP harmonic matrices into LTI HSS formulation
    
    n       = DSA.info.n;
    m       = DSA.info.m;
    p       = DSA.info.p;
    d       = DSA.info.d;
    
    h_t     = DSA.FLIFT.h_t; % maximum truncation order
    nh_t    = 2*h_t+1;       % total number of components
    
    DSA.calc.h_t  = h_t;
    DSA.calc.nh_t = nh_t;
    
    % -------------------------------------------------------------------------
    % definition of lifted names of states, to be used during participation factor
    % analysis:
    
    varnames = DSA.vars.smsi_states;
    smsi_states_lifted = cell(n*nh_t, 1);
    idx = 0;
    for h = -h_t:h_t
        for vi = 1:n
            idx = idx + 1;
            vname = varnames{vi};
            smsi_states_lifted{idx} = sprintf('%s{%+d}', vname, h);
        end
    end
    
    DSA.vars.names.smsi_states_lifted = smsi_states_lifted;
    
    % -------------------------------------------------------------------------
    % From 3D harmonic matrices to 2D lifted matrices
    
    if ~DSA.data.delay.is_delayed
        [TA, TB, TC, TD, TN] = DSA_LTP2HSS(DSA);
        TA_min_TN = TA - TN;
        HSS_SYS = ss(TA_min_TN, TB, TC, TD);

        LTI.mat.A = TA_min_TN;
        LTI.mat.B = TB;
        LTI.mat.C = TC;
        LTI.mat.D = TD;

    else
        [TA, TB1, TB2, TC1, TC2, TD11, TD12, TD21, TD22, TT_d, TM, TN] = DSA_LTP2DHSS(DSA);
        TA_min_TN = TA - TN;
        
        if DSA.FLIFT.cancel_delays % in case you want to easily set all delays to zero.
            TA_min_TN   = TA_min_TN + TB2  * ((eye(d*nh_t)-TD22)\TC2); 
            TB          = TB1       + TB2  * ((eye(d*nh_t)-TD22)\TD21);
            TC          = TC1       + TD12 * ((eye(d*nh_t)-TD22)\TC2);
            TD          = TD11      + TD12 * ((eye(d*nh_t)-TD22)\TD21);
            
            HSS_SYS = ss(TA_min_TN, TB, TC, TD);

            LTI.mat.A = TA_min_TN;
            LTI.mat.B = TB;
            LTI.mat.C = TC;
            LTI.mat.D = TD;
        else
            HSS_SYS = setDelayModel(TA_min_TN, TB1, TB2*TM, TC1, TC2, TD11, TD12*TM, TD21, TD22*TM, TT_d);

            LTI.mat.A   = TA_min_TN;
            LTI.mat.B1  = TB1;
            LTI.mat.B2  = TB2*TM;
            LTI.mat.C1  = TC1;
            LTI.mat.C2  = TC2;
            LTI.mat.D11 = TD11;
            LTI.mat.D12 = TD12*TM;
            LTI.mat.D21 = TD21;
            LTI.mat.D22 = TD22*TM;
            LTI.mat.T_d = TT_d;
        end
    end
    
    DSA.LTI = LTI;
    DSA.LTI.HSS_SYS = HSS_SYS;
    
    % display information about the system :
    disp('--- DSA_get_HSS: the HSS model is defined as DSA.LTI.HSS_SYS')
    fprintf('    ')
    size(HSS_SYS)
    fprintf('    p*nh_t = %2d x %2d = %3d outputs\n', p, nh_t, p*nh_t)
    fprintf('    m*nh_t = %2d x %2d = %3d inputs\n' , m, nh_t, m*nh_t)
    fprintf('    n*nh_t = %2d x %2d = %3d states\n' , n, nh_t, n*nh_t)
    if hasInternalDelay(HSS_SYS)
        disp('    The system has internal (pure) delay(s)')
    else
        disp('    The system does not have internal (pure) delay(s)')
    end
end
