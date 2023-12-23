function x = DSA_HD2TD(X, t, f_1)
%% x = DSA_HD2TD(X, t, f_1)
%
% This function implements the numerical evaluation of a Fourier series.
% More specifically, it transforms Fourier coefficients into a time-domain
% waveform based on the provided time vector and fundamental frequency:
%     t must be a column vector or a scalar [s]
%     X must be a column vector of Fourier coefficients (from - to + freqs)
%     f_1 is the fundamental frequency [Hz]
% 
% Author: Philippe De Rua
% Updated: 2023-06-06

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

if size(t, 2) > 1 % then t is a row vector
    % make it a column vector
    t = t.';
end
if size(X, 2) > 1 % then X is a row vector
    % make it a column vector
    X = X.';
end

x = real(exp(1j*2*pi*f_1 * t*(-((length(X)-1)/2):((length(X)-1)/2))) * X);
end