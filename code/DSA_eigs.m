function DSA = DSA_eigs(DSA)
%% DSA = DSA_eigs(DSA)
%
% Dispatch function DSA_eigs runs eigenvalues calculations for non-delayed
% systems and characteristic roots calculations for delayed systems. In the
% latter case, only the Arnoldi Algorithm is interfaced with the present
% toolbox. However, please be aware of the fact that there exist other
% routines which you may want to try as well.
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

if DSA.data.delay.is_delayed
    % calculation of DDE characteristic roots

    str = ['-- DSA_eigs: started calculation of DDE characteristic roots\n', ...
    '\n', ...
    '   This calculation relies on a separate toolbox developed\n', ...
    '   and made available by the KU Leuven Computer Science Department.\n', ...
    '   Software webpage:\n', ...
    '           https://twr.cs.kuleuven.be/research/software/delay-control/\n', ...
    '           Accessed: 2021\n', ...
    '   Contact Person:\n', ...
    '           Wim Michiels\n', ...
    '           https://people.cs.kuleuven.be/~wim.michiels/\n'];
    fprintf(str);

    DDE_start = tic;
    DSA = DSA_DDE_chroots(DSA);
    DDE_end = toc(DDE_start);
    
    disp('Total elapsed time: ')
    disp(datestr(datenum(0,0,0,0,0,DDE_end), 'HH:MM:SS'))
else
    % calculation of ODE eigenvalues
    % ddt(x(t)) = A x(t)
    
    disp('-- DSA_eigs: started calculation of ODE eigenvalues')
    
    [V, D, W] = eig(DSA.LTI.mat.A); % for HSS system, LTI.mat.A = TA-TN
    
    EIGS = diag(D);
    DSA.LTI.EIGS = EIGS;
    DSA.LTI.RVEC = V; % columns are the right eigenvectors
    DSA.LTI.LVEC = W; % *columns* (not rows) are the left eigenvectors
    
    disp('-- DSA_eigs: done')
end

%% figure

EIGS = DSA.LTI.EIGS;

figure
hold on
grid on
plot(real(EIGS), imag(EIGS)/(2*pi), '.', 'markersize', 12)
if strcmp(get(groot, 'defaultTextInterpreter'), 'latex')
    xlabel('real part: $\sigma$ [rad/s]', 'fontsize', 14)
    ylabel('imag part: $\omega/(2\pi)$ [Hz]', 'fontsize', 14)
else
    xlabel('real part: \sigma [rad/s]', 'fontsize', 14)
    ylabel('imag part: \omega/(2\pi) [Hz]', 'fontsize', 14)
end
yl = ylim;
if abs(yl(1)-yl(2)) < 5e-4
   ylim(yl(1) + [-0.5 +0.5])
end
xl = xlim;
if abs(xl(1)-xl(2)) < 5e-4
   xlim(xl(1) + [-0.5 +0.5])
end
sgtitle('eigenvalues/characteristic roots')

end