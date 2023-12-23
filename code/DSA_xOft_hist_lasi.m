function x_hist_lasi = DSA_xOft_hist_lasi(t, varnames, HD_struct)
%% x_hist_lasi = DSA_xOft_hist_lasi(t, varnames, HD_struct)
%
% Gives the time-domain waveform history of the variables in varnames,
% based on the Fourier coefficients stored in the HD structure HD_struct.
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

    sz = size(varnames, 1);
    x_hist_lasi = zeros(sz, length(t));
    for vi = 1:sz
        x_hist_lasi(vi, :) = DSA_HD2TD(HD_struct.(varnames{vi}), t, HD_struct.f_b); % in time domain
    end
end