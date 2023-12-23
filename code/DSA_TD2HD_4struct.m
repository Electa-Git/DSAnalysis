function [HD_out, FullHD_out] = DSA_TD2HD_4struct(TD_in, varnames, f_s, h_max)
%% [HD_out, FullHD_out] = DSA_TD2HD_4struct(TD_in, varnames, f_s, h_max)
% 
% Implementation of the DFT - Discrete Fourier Transform
% The function takes time-domain waveform x_TD as input. Its corresponding
% sampling frequency must be specified as f_s. The maximum harmonic index
% up to which the Fourier coefficients are of interest (i.e. the requested
% maximum harmonic rank) is given in h_max.
%
% This function relies on DSA_DFT and automates the calculation for
% structures.
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

% extracting the frequency vector and max availble harmonic :
[~, f_full, h_max_available] = DSA_DFT(TD_in.(varnames{1}), f_s, false, true, '');

% automatised check and mappings :
if isnan(h_max)
    h_max = h_max_available;
elseif h_max > h_max_available
   warning('The requested h_max is larger than the maximum available harmonic index. Try increasing the sampling frequency. In the meantime, extraction of harmonic data is done up to the maximum available harmonic index.') 
   h_max = h_max_available;
end

zeroHz_index = find(f_full == 0);
if isempty(zeroHz_index)
    zeroHz_index = find(abs(f_full) == min(abs(f_full)));
end

h2i_full = @(h) h + zeroHz_index;
freq_map = h2i_full(-h_max:h_max);

FullHD_out.f   = f_full;
FullHD_out.f_b = f_full(h2i_full(1));

HD_out.f   = f_full(freq_map);
HD_out.f_b = f_full(h2i_full(1));

% loop on every variable :
for k = 1:size(varnames)
    % getting the full spectrum of the variable by DISCRETE FOURIER TRANSFORM :
    [harmonic_vector, ~, ~] = DSA_DFT(TD_in.(varnames{k}), f_s);
    FullHD_out.(varnames{k}) = harmonic_vector;

    % extracting only the harmonics over -h_max:h_max from the full spectrum :
    HD_out.(varnames{k}) = FullHD_out.(varnames{k})(freq_map);
end
end