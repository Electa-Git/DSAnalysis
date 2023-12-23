function DSA = DSA_new(sstm, inst)
%% DSA = DSA_new(sstm, inst)
% Creation of a new DSA structure.
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

% assigning system and instance information:
DSA.sstm = sstm;
DSA.inst = inst;

% creating empty structures for storage of data and results:
DSA.info = struct();
DSA.fcts = struct();
DSA.calc = struct();
DSA.TD   = struct();
DSA.HD   = struct();

% --------------------------------------------------------------------
% default values for system's first initialisation:

% h_m             : maximum harmonic rank for periodic trajectory calculations
% jacobian_type   : 'numerical' or 'analytical'
% jacobian_reset  : 1 to re-generate the analytical jacobian via symbolic 
%                   calculation. This parameter is ignored if jacobian_type
%                   = 'numerical'.
% jacobian_eigs   : 1 to enable the calculation of jacobian eigenvalues.
%                   This parameter is ignored for systems with pure delays.
% def_inputs_from : input variables defined from 'user', 'source', or 'solution'
% init_guess_from : initial guess defined from 'user', 'source', or 'solution'
% HD_source       : HD source name for initialisation, usually 'HD_COLL_1P'
% algo_tol        : algorithm tolerance
% max_iter        : maximum number of iterations
% L_DISP          : number of sample points per fundamental period for the
%                   display of the waveforms

DSA.COLL.h_m             = 6;
DSA.COLL.jacobian_type   = 'analytical';
DSA.COLL.jacobian_reset  = 1;
DSA.COLL.jacobian_eigs   = 1;
DSA.COLL.def_inputs_from = 'user';
DSA.COLL.init_guess_from = 'solution';
DSA.COLL.HD_source       = 'HD_COLL_1P';
DSA.COLL.algo_tol        = 1e-9;
DSA.COLL.max_iter        = inf;
DSA.COLL.L_DISP          = 1625;

% XP        : Number of fundamental periods for Fourier analysis
% NP        : Number of fundamental periods for numerical integration
% init_from : 'solution' or 'source' or 'user' or 'zero'
% HD_source : HD source for initialisation
% sim_tol   : integration tolerance
% solver    : one of matlab's differential equation solvers. The toolbox
%             automatically uses dde23 for delayed systems.

DSA.TDSIM.XP        = 2;
DSA.TDSIM.tspan     = [0, 0.062];
DSA.TDSIM.T_s       = 10e-6;
DSA.TDSIM.init_from = 'source';
DSA.TDSIM.HD_source = 'HD_COLL_1P';
DSA.TDSIM.sim_tol   = 1e-8;
DSA.TDSIM.solver    = 'ode15s';

% h_t           : maximum harmonic index for the frequency lifting method
% h_f           : forced periodic rank, generally kept equal to h_m. Leave
%                 empty as [] to automatically set h_f equal to h_m.
% cancel_delays : set to true to rewrite pure delays to zero.

DSA.FLIFT.h_t = 3;
DSA.FLIFT.h_f = [];
DSA.FLIFT.cancel_delays = 0;

% --------------------------------------------------------------------
% defining function folders and adding their paths

path_sstm = ['.\sstms\' sstm];
path_inst = ['.\sstms\' sstm '\' inst];

addpath(path_sstm)
addpath(path_inst)
addpath('.\sstms')
addpath(genpath('.\code')) % adds path *and* subfolders

DSA.info.path_sstm = path_sstm;
DSA.info.path_inst = path_inst;

% --------------------------------------------------------------------
% defining generic functions

str_DATA = [path_sstm '\' sstm          '_DATA.m'];
str_VARS = [path_inst '\' sstm '_' inst '_VARS.m'];
str_HD   = [path_inst '\' sstm '_' inst '_HD.m'];
str_uOft = [path_inst '\' sstm '_' inst '_uOft.m'];
str_EQNS = [path_inst '\' sstm '_' inst '_EQNS.m'];

DSA.fcts.DATA         = eval(['@(DSA)     ' sstm          '_DATA  (DSA)']);      % function
DSA.fcts.VARS         = eval(['@(DSA)     ' sstm '_' inst '_VARS  (DSA)']);      % function
DSA.fcts.HD           = eval(['@(DSA, HDi)' sstm '_' inst '_HD    (DSA, HDi)']); % function
DSA.fcts.uOft         = eval(['@(t, DSA)  ' sstm '_' inst '_uOft  (t, DSA)']);   % function
DSA.fcts.EQNS         = [sstm '_' inst '_EQNS'];         % script
DSA.fcts.COLL_JAC     = [sstm '_' inst '_COLL_JAC'];     % script
DSA.fcts.COLL_JAC_DLY = [sstm '_' inst '_COLL_JAC_DLY']; % script
DSA.fcts.open         = @(varargin) DSA_open_project(str_DATA, str_VARS, str_HD, str_uOft, str_EQNS, varargin);
DSA.fcts.save         = @(varargin) DSA_save_project(path_inst, varargin);

% --------------------------------------------------------------------
% populating case-specific fields and getting symbolic nonlinear system

DSA = DSA_set(DSA); % set numerical data, variables lists and calculation parameters

disp('A new DSA structure has been created:')
disp(['        System           : ' sstm])
disp(['        Instance         : ' inst])
disp(['        State variables  : ' num2str(DSA.info.N)])
disp(['        Input variables  : ' num2str(DSA.info.M)])
disp(['        Output variables : ' num2str(DSA.info.P)])
disp(['        Delay variables  : ' num2str(DSA.info.D)])
if (DSA.info.D>0) && ~DSA.data.delay.is_delayed
    disp('        (all delays are zero)')
end
end


