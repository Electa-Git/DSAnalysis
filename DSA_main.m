%% DSA: Dynamic System Analysis
% 
% Getting started:
% 
% (1) choose system [sstm] and instance [inst]
% (2) generate a new DSA structure with DSA_new(sstm, inst)

sstm = 'RLC';
inst = 'Ch3';

% sstm = 'MMC';
% inst = 'CCC_1PH';
% inst = 'CCC_1PH_DLY';
% inst = 'CCC_3PH';
% inst = 'CCC_3PH_DLY';

% sstm = 'NET';
% inst = 'GFM_MMC_OWF_3TLC_DLY';
% inst = 'HVDC_GFL_MMC_DLY';
% inst = 'HVDC_GFM_MMC_OWF_3TLC_DLY';
% inst = 'OWF_1TLC_DLY';
% inst = 'OWF_3TLC_DLY';

% getting new DSA structure with default values:
DSA = DSA_new(sstm, inst);

%% Set/Open/Save project

% DSA = DSA_set(DSA); % can be called at any point, for instance to update
% numerical parameter values. COLL and TDSIM always call DSA_set as first
% step before starting their calculations, thereby making sure parameter
% values are up-to-date. DSA_set is also called by DSA_new.

% DSA.fcts.open()
% DSA.fcts.open('EQNS')
% DSA.fcts.open('DATA')
% DSA.fcts.open('HD')
% DSA.fcts.open('uOft')

% DSA.fcts.save()
% DSA.fcts.save('example')

%% Find system steady-state periodic trajectory with COLL
% COLL: Fourier-based collocation method

DSA.COLL.h_m             = 1;
DSA.COLL.jacobian_type   = 'numerical'; % 'numerical' or 'analytical'
DSA.COLL.jacobian_reset  = 0;
DSA.COLL.def_inputs_from = 'user'; % 'user', 'source' or 'solution'
DSA.COLL.init_guess_from = 'user'; % 'user', 'source' or 'solution'
DSA.COLL.HD_source       = 'HD_COLL_1P';
DSA.COLL.algo_tol        = 1e-9;
DSA.COLL.L_DISP          = 1625;

DSA = DSA_run_COLL(DSA);

% DSA_plot(DSA)

%% Direct time-domain integration (nonlinear system)
% Simulation of the nonlinear system

DSA.TDSIM.T_s       = 10e-6;
DSA.TDSIM.XP        = 3;
DSA.TDSIM.tspan     = [0, 0.02];
DSA.TDSIM.init_from = 'source'; % 'user', 'source', 'zero' or 'solution'
DSA.TDSIM.HD_source = 'HD_COLL_1P';
DSA.TDSIM.abs_tol   = 1e-6;
DSA.TDSIM.rel_tol   = 1e-6;
DSA.TDSIM.solver    = 'ode15s'; % e.g. 'ode45', 'ode15s', 'ode23tb',...

DSA = DSA_run_TDSIM(DSA);

% DSA_plot(DSA)

%% positive-negative-zero sequences calculations
% only for three-phase systems (abc variables)

HD_source = 'HD_COLL_1P';
TD_source = 'TD_COLL_1P';

% example:
vnames = {'v_g';
          'i_s';
          'v_c_u';
          'v_c_l';
         {'v_c', 'ref'}};
DSA = DSA_sequences(DSA, HD_source, TD_source, vnames);

%% Symbolic manipulations

% get symbolic nonlinear system:
% (note that the symbolic nonlinear system may have been calculated when
% running the collocation method, in case an analytical jacobian is used
% AND has been reset)
DSA = DSA_get_NLSS(DSA);

% symbolic linearisation of the nonlinear system, to obtain symbolic
% expressions of the state-space matrix coefficients:
DSA = DSA_linearise(DSA);

%% Fourier analysis of periodic coefficients of LTP system

DSA.LTP.HD_source = 'HD_COLL_1P';
DSA = DSA_get_3D_LTP(DSA);

%% Applying frequency-lifting

DSA.FLIFT.h_t = 1;
DSA.FLIFT.h_f = DSA.COLL.h_m;
DSA = DSA_get_HSS(DSA);

%% Frequency response calculation

% TFid        : transfer function identifier
% ivarname    : input smsi variable (available in DSA.vars.smsi_input)
% ovarname    : name of output smsi variable
% sign_factor : +1 or -1 to change the sign convention
% h_in        : input harmonic channel
% h_shift     : harmonic shift
% f_values    : vector of frequency values

TFid = 'SISO_example';
DSA.FRESP.(TFid).ivarname     = 'v_d';
DSA.FRESP.(TFid).ovarname     = 'i_d';
DSA.FRESP.(TFid).sign_factor  = 1;
DSA.FRESP.(TFid).h_in         = 0;
DSA.FRESP.(TFid).h_shift      = 0;
DSA.FRESP.(TFid).f_values     = logspace(-2, 4, 1001);
DSA = DSA_run_FRESP(DSA, TFid);

TFid = 'MIMO_3PH_example';
DSA.FRESP.(TFid).ivarname     = 'v_g';
DSA.FRESP.(TFid).ovarname     = 'i_s';
DSA.FRESP.(TFid).sign_factor  = -1;
DSA.FRESP.(TFid).h_in         = 0;
DSA.FRESP.(TFid).h_shift      = 0;
DSA.FRESP.(TFid).f_values     = logspace(-2, 4, 1001);
DSA = DSA_run_FRESP_3PH(DSA, TFid);

%% Analysis eigenvalues

% options related to DDE characteristic roots (disregarded if the system
% does not have delays):
% nb_eig_fctr      : nb of roots to be extracted, as a fraction of the 
%                    total nb of states (e.g. 0.8, 1.0, 1.2,...)
% convergence_test : set to 1 to check which roots have converged
% convergence_eps  : error defining root convergence, e.g. 1e-10

DSA.LTI.DDE.nb_eig_fctr      = 1;
DSA.LTI.DDE.convergence_test = 1;
DSA.LTI.DDE.convergence_eps  = 1e-10;

DSA = DSA_eigs(DSA);

%% Participation factors analysis
% Function DSA_partfactors uses its own pfa data structure.

% EIGS     : eigenvalues
% RVEC     : right eigenvectors
% LVEC     : left eigenvectors
% VARS     : cell of states names
% MODES    : requested modes (first, run with [])
% XMIN_TXT : indices are displayed for eigenvalues of real part larger than
%            XMIN_TXT
% YMIN_TXT : indices are displayed for eigenvalues of imaginary part larger
%            than YMIN_TXT

pfa.EIGS     = DSA.LTI.EIGS;
pfa.RVEC     = DSA.LTI.RVEC;
pfa.LVEC     = DSA.LTI.LVEC;
pfa.VARS     = DSA.vars.names.smsi_states_lifted;
pfa.MODES    = [];
pfa.XMIN_TXT = -inf;
pfa.YMIN_TXT = -inf;

DSA_partfactors(pfa);

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