function DSA_plot(DSA)
%% DSA_plot(DSA)
%
% Function DSA_plot helps you quickly plot time-domain waveforms and
% harmonic-domain spectra from the COLL (collocation) method and the TDSIM
% (time-domain numerical simulation) method.
%
% Author: Philippe De Rua
% Updated: 2023-11-28

%#ok<*NASGU> 

%% (1) run this code snippet to get the list of all available variables
%  (2) copy-paste the list to (re)define the 'sets' variable below
%  (3) choose your TD or HD source
%  (4) run the code to plot the variables
    
DSA_dispvars(DSA)

%% Plot TD waveforms and HD spectra

% TD_sources cell content: source name, line type, legend text, density reduction
TD_sources = {
'TD_COLL_1P',   {'-'},   'COLL 1P',  1;
% 'TD_TDSIM_1P',  {'-'},  'TDSIM 1P', 1;
% 'TD_TDSIM_XP',  {'-'},  'TDSIM XP', 1;
% 'TD_TDSIM_NP',  {'-'},  'TDSIM NP', 1;
}; 

% HD_sources cell content: source name, legend text
HD_sources = {
% 'HD_COLL_1P',  'COLL 1P';
% 'HD_TDSIM_1P', 'TDSIM 1P';
% 'HD_TDSIM_XP', 'TDSIM XP';
% 'FHD_TDSIM_XP', 'TDSIM XP';
};

sets = {{
    % to be filled in.
}};

run('DSA_run_plot')

end
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