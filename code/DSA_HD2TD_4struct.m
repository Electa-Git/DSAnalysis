function TD_out = DSA_HD2TD_4struct(HD_in, fieldnames, t)
%% TD_out = DSA_HD2TD_4struct(HD_in, fieldnames, t)
%
% Function DSA_HD2TD_4struct transforms harmonic vectors into TD waveforms
% waveform based on the provided time vector and fundamental frequency. It
% is generalised to directly process whole structures.
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

f_b = HD_in.f_b;
TD_out.t = t;

% loop on every variable:
for k = 1:size(fieldnames, 1)
    TD_out.(fieldnames{k, 1}) = DSA_HD2TD(HD_in.(fieldnames{k, 1}), t, f_b);
end
end