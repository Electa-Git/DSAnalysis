function output = DSA_dxdt_lasi(t, x, DSA)
%% output = DSA_dxdt_lasi(t, x, DSA)
%
% This function returns the dxdt vector of ODE for their numerical
% integration.
% This function is also used to retrieve output values over complete time
% spans once the numerical integration is finished.
%
% Author: Philippe De Rua
% Updated: 2023-06-09

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

% fetching input variables from uOft:
inputs = DSA.fcts.uOft(t, DSA);

% fetching state variables in vector x from their index:
if DSA.calc.extract_all
    x_source = @(idx) x(:, idx);
else
    x_source = x;
end
for idx = 1:DSA.info.N
    states.(DSA.vars.lasi_states{idx}) = x_source(idx); %#ok<STRNU>
end

% retrieving dxdt:
get_values  = 1;
run_algebra = 1;
run_dxdt    = ~DSA.calc.extract_all;
is_delayed  = false;
param       = DSA.data.param;

if DSA.info.D > 0
    % delayed variables are requested by the EQNS file, even though they
    % are not used. They are thus defined as empty vectors in structure
    % 'delayed':
    for idx = 1:DSA.info.D
        delayed.(DSA.vars.lasi_delayed{idx}) = []; %#ok<STRNU>
    end
end

run(DSA.fcts.EQNS)

% extracting output variables
if DSA.calc.extract_all
    varnames = [DSA.vars.lasi_states;
                DSA.vars.lasi_inputs;
                DSA.vars.lasi_outputs];
    for vi = 1:size(varnames)
        vname = varnames{vi};
        TD_DE.(vname) = eval(vname);
    end
    
    if DSA.info.D > 0
        % If the system allows for pure delays but all delays are set to
        % zero, then the delayed variables are not defined.
        % The solution to this problem is to set them equal to
        % their todelay counterparts:
        for vi = 1:DSA.info.D
            vname_todelay = DSA.vars.lasi_todelay{vi};
            vname_delayed = DSA.vars.lasi_delayed{vi};
            TD_DE.(vname_todelay) = eval(vname_todelay);
            TD_DE.(vname_delayed) = TD_DE.(vname_todelay);
        end
    end
    
    output.TD_DE = TD_DE;
else
    output = dxdt;
end
end