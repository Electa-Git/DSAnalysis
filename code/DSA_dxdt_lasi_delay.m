function output = DSA_dxdt_lasi_delay(t, x, X_dly, DSA)
%% output = DSA_dxdt_lasi_delay(t, x, X_dly, DSA)
%
% This function returns the dxdt vector of DDEs for dde23.
% The presence of multiple time delays is handled automatically.
% This function is also used to retrieve output values over complete time
% spans once the numerical integration is finished.
%
% When integrating DDEs via dde23:
% x is a column vector. This column corresponds to a time t, rows
% correspond to the state variables.
% X_dly is a matrix. Columns correspond to t-T_d1, t-T_d2, etc. Rows
% correspond to the state-variables.
% 
% When extracting data, x is a matrix. The time dimension is along the 
% second dimension, so the columns correspond to different time steps.
% X_dly is a 3D matrix. Columns are still according to the various values
% of delays, while the time dimension is now along the third dimension.
%
% Author: Philippe De Rua
% Updated: 2023-06-07

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

delay = DSA.data.delay;
param = DSA.data.param;
uOft  = @(t) DSA.fcts.uOft(t, DSA);

% ---------------------------------------------------------------------
% handling delays:

% The calculation of the states derivatives at time t involves quantities
% that depend on delayed variables, i.e. on states and inputs as they were
% at previous time steps, e.g. at time t-Td.
% For this reason, this dxdt function first enters a loop where variables
% can be calculated as they were at previous time steps.
% While delayed states are provided by dde23 in X_dly, the delayed inputs
% must be retrieved from uOft manually.
% At the end of each passage in this loop, delayed variables are stored.

% delayed variables are being calculated but are not known yet:
for idx = 1:DSA.info.D
    delayed.(DSA.vars.lasi_delayed{idx}) = [];
end

% determining the necessary calculations, where delays equal to zero are
% evaluated last:
if any(delay.unique_T_d == 0)
    needed_calc = [1:length(delay.unique_nonzero_T_d), 0];
else
    needed_calc = 1:length(delay.unique_nonzero_T_d);
end

get_values  = 1;
run_algebra = 1;
run_dxdt    = 0;
is_delayed  = 0; % to calculate the delayed variables at the tested times.

for calc = needed_calc
    if calc > 0
        % calculation of delayed signals
        this_delay = delay.unique_nonzero_T_d(calc);
        x_source   = permute(X_dly(:, calc, :), [1, 3, 2]); % "permute" makes a 2D matrix out of it.
        inputs     = uOft(t - this_delay);
    else
        % calculation of non-delayed signals
        this_delay = 0;
        x_source   = x;
        inputs     = uOft(t);
    end
    
    % fetching state variables from their index:
    if size(x_source, 2) ~= 1
        for idx = 1:DSA.info.N
            states.(DSA.vars.lasi_states{idx}) = x_source(idx,:).';
        end
    else
        for idx = 1:DSA.info.N
            states.(DSA.vars.lasi_states{idx}) = x_source(idx);
        end
    end

    run(DSA.fcts.EQNS)
    
    % ---------------------------------------------------------------------
    % storing delayed versions of to-delay variables:
    % ---------------------------------------------------------------------
    % example:
    % the code just ran with Td = 10 ms. This is the delay of variable v.
    % Following this run, v has a "delayed" value, which is stored:
    % delayed.v = v;
    
    for idx = 1:DSA.info.D
        lasi_todelay_varname = DSA.vars.lasi_todelay{idx};
        if delay.(lasi_todelay_varname) == this_delay
            lasi_delayed_varname = DSA.vars.lasi_delayed{idx};
            delayed.(lasi_delayed_varname) = eval(lasi_todelay_varname);
        end
    end
end

%% Now that all delayed variables have been stored, they are available via
% the delayed structure.
% These delayed values can now be used for all calculations involving
% delayed variables.

% -------------------------------------------------------------------------
% calculations using delayed signals:
% -------------------------------------------------------------------------

if ~any(delay.unique_T_d == 0)
    % in this case, values at time t are not yet available in inputs and
    % states structures.
    
    % fetching input variables:
    inputs = uOft(t);
    
    if DSA.calc.extract_all
        x_source = @(idx) x(idx, :).';
    else
        x_source = x;
    end

    % fetching state variables from their index:
    for idx = 1:DSA.info.N
        states.(DSA.vars.lasi_states{idx}) = x_source(idx);
    end
end

get_values  = 1;
run_algebra = 1;
run_dxdt    = ~DSA.calc.extract_all;
is_delayed  = 1; % this time, working with the delayed variables.

run(DSA.fcts.EQNS)

if DSA.calc.extract_all
    varnames = DSA.vars.lasi_all;
    for vi = 1:length(varnames)
        vname = varnames{vi};
        TD_DE.(vname) = eval(vname);
    end

    output.TD_DE = TD_DE;
else
    output = dxdt;
end
% -------------------------------------------------------------------------
end