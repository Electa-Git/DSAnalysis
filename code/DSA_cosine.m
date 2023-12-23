function output = DSA_cosine(Mag, Phase)
%% output = DSA_cosine(Mag, Phase)
%
% Function DSA_cosine defines the Fourier coefficients of a cosine wave
% based on the amplitude (any units) and phase angle (in degrees!).
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

output = [Mag*exp(-1j*Phase/180*pi)/2 Mag*exp(1j*Phase/180*pi)/2];

end