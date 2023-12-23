function [x_HD, f, f_b] = DSA_TD2HD(x_TD, f_s, h_max)
%% [x_HD, f, f_b] = DSA_TD2HD(x_TD, f_s, h_max)
%
% Implementation of the DFT - Discrete Fourier Transform
% The function takes time-domain waveform x_TD as input. Its corresponding
% sampling frequency must be specified as f_s. The maximum harmonic index
% up to which the Fourier coefficients are of interest (i.e. the requested
% maximum harmonic rank) is given in h_max.
%
% The function returns the vector of Fourier coefficients in x_HD, as well
% as the vector of corresponding frequencies in f and the value of the base
% frequency f_b.
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

% getting the full spectrum of the variable by DISCRETE FOURIER TRANSFORM :
[x_fullHD, f_full, h_max_available] = DSA_DFT(x_TD, f_s, 0, 0, '');

% automatised check and mappings :
if isnan(h_max)
    h_max = h_max_available;
elseif h_max > h_max_available
   error('PROBLEM: requested h_max is too large. The sampling frequency must be increased.') 
end

DC_freq_index = find(f_full == 0);

if isempty(DC_freq_index)
    error('PROBLEM: no 0 Hz frequency found.') 
end

h2i_full = @(h) h + DC_freq_index;
freq_map = h2i_full(-h_max:h_max);

% extracting only the harmonics over -h_max:h_max from the full spectrum :
x_HD    = x_fullHD(freq_map).';
f       = f_full(freq_map);
f_b     = f_full(h2i_full(1));
end