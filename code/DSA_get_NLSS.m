function DSA = DSA_get_NLSS(DSA)
%% DSA = DSA_get_NLSS(DSA)
%
% This function creates the symbolic nonlinear system from file EQNS.
% 
% Author: Philippe De Rua
% Updated: 2023-06-08

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

%% Definition of the variables

% -------------------------------------------------------------------------
% Making all variables and parameters symbolic:

to_symbolic = [
    DSA.vars.lasi_inputs;
    DSA.vars.lasi_states;
    DSA.vars.lasi_delayed;
    fields(DSA.data.param)];

for idx = 1:size(to_symbolic,1)
    eval([to_symbolic{idx} ' = sym(''' to_symbolic{idx} ''');']);
end

%% Symbolic calculations

%#ok<*NASGU>
%#ok<*NODEF>

% -------------------------------------------------------------------------
% dxdt :
% -------------------------------------------------------------------------

get_values  = 0; % no need for numerical values in symbolic calculations
run_algebra = 1; % --> running algebraic relationships leading to f(x,u)
run_dxdt    = 1; % --> getting the expression of dxdt = f(x,u)
is_delayed  = DSA.data.delay.is_delayed;

run(DSA.fcts.EQNS)

% -------------------------------------------------------------------------
% outputs expressions:
% -------------------------------------------------------------------------

varnames = DSA.vars.smsi_outputs;
if isempty(varnames)
    warning('DSA_get_NLSS: smsi_outputs is empty.')
    y_smsi = [];
else
    y_smsi = eval(cell2sym(varnames));
end

% -------------------------------------------------------------------------
% to-delay variable expressions:
% -------------------------------------------------------------------------

varnames = DSA.vars.lasi_todelay;
if isempty(varnames)
    z_lasi = [];
else
    z_lasi = eval(cell2sym(varnames));
end

varnames = DSA.vars.smsi_todelay;
if isempty(varnames)
    z_smsi = [];
else
    z_smsi = eval(cell2sym(varnames));
end

% -------------------------------------------------------------------------

NLSS.dxdt   = dxdt;
NLSS.lasi.z = z_lasi;
NLSS.smsi.y = y_smsi;
NLSS.smsi.z = z_smsi;
DSA.NLSS    = NLSS;

disp('--- DSA_get_NLSS: symbolic nonlinear system defined in DSA.NLSS')
end
