function [Y, f, h_max_available] = DSA_DFT(varargin)
%% [Y, f, h_max_available] = DSA_DFT(varargin)
% [Y, f, h_max_available] = DSA_DFT(y, f_s)
% [Y, f, h_max_available] = DSA_DFT(y, f_s, plot_figs, verbose, name_str)
% 
% DSA_DFT essentially calculates the Fourier coefficients of a given
% time-domain waveform. It relies on the fft. It automatically handles the
% calculation of the base frequency, which is determined from the number of
% elements in time-domain vector y and the sampling frequency f_s.
%
% INPUTS:
% y is the time-domain signal;
% f_s is its sampling frequency;
% plot_figs specifies whether to plot the figures or not.
% name_str is used as a title to plot the periodic signal.
% 
% OUTPUTS:
% Y is the vector of Fourier coefficients of y, from - to + frequencies;
% f is the zero-centered vector of frequencies corresponding to y.
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

switch nargin
    case 5
        y           = varargin{1};
        f_s         = varargin{2};
        plot_figs   = varargin{3};
        verbose     = varargin{4};
        name_str    = varargin{5};
    case 2
        y           = varargin{1};
        f_s         = varargin{2};
        plot_figs   = false;
        verbose     = false;
        name_str    = '';
    otherwise
        error('incorrect number of arguments')
end

L = length(y);
Y = fftshift(fft(y)/L);
df = f_s/L; % valid for both even and odd L
f = (-f_s/2:df:f_s/2-df) + mod(L,2)*df/2; % valid for both even and odd L
if mod(L,2) == 0 % even length : non-symmetrical harmonic content around 0
    h_max_available = (L-2) / 2;
else % odd length : symmetrical harmonic content around 0
    h_max_available = (L-1) / 2;
end

if verbose
    disp(' ')
    disp('--------------------')
    disp('  Fourier Analysis  ')
    disp('--------------------')
    disp(['L           = ' num2str(L) ' samples'])
    disp(['f_s         = ' num2str(f_s) ' Hz'])
    disp(['df          = ' num2str(df) ' Hz'])
    disp(['f_min       = ' num2str(f(1)) ' Hz'])
    disp(['f_max       = ' num2str(f(end)) ' Hz'])
    % disp(['h_max_avail = ' num2str(h_max_available)])
    disp(' ')
end

if plot_figs
    figure;
    ax1 = subplot(2, 1, 1);
        stem(f, abs(Y), '.')
        xlim(f([1, end]))
        if max(abs(Y)) > 0
            ylim([0, max(abs(Y))*1.1])
        end
        xlabel('frequency [Hz]')
        ylabel('amplitude [unit]')
        grid
        title(name_str)
    ax2 = subplot(2, 1, 2);
        stem(f, angle(Y)*180/pi, '.')
        xlim(f([1, end]))
        ylim([-190, +190])
        xlabel('frequency [Hz]')
        ylabel('phase [deg]')
        grid
    linkaxes([ax1, ax2], 'x')
end
end