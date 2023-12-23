function DSA = DSA_linearise(DSA)
%% DSA = DSA_linearise(DSA)
%
% This function linearises the nonlinear system, i.e. calculates the 
% symbolic matrix coefficients of the corresponding linear state-space
% representation.
%
% The linearisation is carried out between smsi_inputs and smsi_outputs,
% which are potentially different from lasi_inputs and lasi_outputs.
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

disp('--- DSA_linearise: symbolic calculation of linear state-space matrix coefficients')

NLSS = DSA.NLSS;
dxdt = NLSS.dxdt;
y    = NLSS.smsi.y;
z    = NLSS.smsi.z;

if ~DSA.data.delay.is_delayed
    local_states = DSA.vars.smsi_states;
    local_inputs = DSA.vars.smsi_inputs;
    
    %                                 derivatives with respect to the :
    mat.A = jacobian(dxdt, cell2sym(local_states)); % STATE variables
    mat.B = jacobian(dxdt, cell2sym(local_inputs)); % INPUT variables
    mat.C = jacobian(   y, cell2sym(local_states)); % STATE variables
    mat.D = jacobian(   y, cell2sym(local_inputs)); % INPUT variables
    
else
    local_states  = DSA.vars.smsi_states;
    local_inputs  = DSA.vars.smsi_inputs;
    local_delayed = DSA.vars.smsi_delayed;
    
    %                                derivatives with respect to the :
    mat.A   = jacobian(dxdt, cell2sym(local_states )); % STATE   variables
    mat.B1  = jacobian(dxdt, cell2sym(local_inputs )); % INPUT   variables
    mat.B2  = jacobian(dxdt, cell2sym(local_delayed)); % DELAYED variables
    mat.C1  = jacobian(   y, cell2sym(local_states )); % STATE   variables
    mat.D11 = jacobian(   y, cell2sym(local_inputs )); % INPUT   variables
    mat.D12 = jacobian(   y, cell2sym(local_delayed)); % DELAYED variables
    mat.C2  = jacobian(   z, cell2sym(local_states )); % STATE   variables
    mat.D21 = jacobian(   z, cell2sym(local_inputs )); % INPUT   variables
    mat.D22 = jacobian(   z, cell2sym(local_delayed)); % DELAYED variables
end

DSA.LNSS.mat = mat;

disp('--- DSA_linearise: symbolic linearised system stored in DSA.LNSS')

end
