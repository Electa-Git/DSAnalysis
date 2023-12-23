function DSA_dispvars(DSA)
%% DSA_dispvars(DSA)
%
% Function DSA_dispvars displays the list of all variables available
% for plot (time-domain waveforms or harmonic-domain spectra.
% See DSA_plot file to learn about its use.
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

fprintf('\nsets = {{\n%% -------------------- %s_%s\n', DSA.sstm, DSA.inst);
local_lasi_all = [DSA.vars.lasi_states;
                  DSA.vars.lasi_inputs;
                  DSA.vars.lasi_outputs;
                  DSA.vars.lasi_todelay;
                  DSA.vars.lasi_delayed];
for k = 1:length(local_lasi_all)
    varname = local_lasi_all{k};
    fprintf('\n%%     ''%s'';', varname);
end
     
fprintf('\n}};\n');
end