function DSA = DSA_DDE_chroots(DSA)
%% DSA = DSA_DDE_chroots(DSA)
%
% Calculation of characteristic roots of homogeneous delay differential
% equations.
% 
%    This calculation relies on a separate toolbox developed
%    and made available by the KU Leuven Computer Science Department.
%    Software webpage:
%            https://twr.cs.kuleuven.be/research/software/delay-control/
%            Accessed: 2021
%    Contact Person:
%            Wim Michiels
%            https://people.cs.kuleuven.be/~wim.michiels/
%
% Author: Philippe De Rua
% Updated: 2023-04-23

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
    
    dly_values = DSA.data.delay.unique_nonzero_T_d;
    mat = DSA.LTI.mat;
    DDE_A   = mat.A;
    DDE_B2  = mat.B2;
    DDE_C2  = mat.C2;
    DDE_D22 = mat.D22;
    DDE_T   = mat.T_d;
    n       = DSA.info.n;
    nh_t    = DSA.calc.nh_t;
    
    if any(any(DDE_D22 ~= 0))
        error('-- DSA_DDE_chroots: the code does not (yet) allow for non-zero D22 matrices')
    end
    
    % The number of state matrices is equal to the number of unique delay
    % values plus one matrix related to delay-free dynamics:
    nb_matrices = length(dly_values) + 1;
    tds_matrices = cell(1, nb_matrices);
    
    % Vector 'to_keep' determines which dynamics are related to which delay
    % values.
    % Calculations related to the no-delay state matrix:
    to_keep = DDE_T == 0;
    tds_matrices{1} = DDE_A + DDE_B2*diag(to_keep)*DDE_C2; 
    % Calculations related to the delay state matrices:
    for iTd = 1:length(dly_values)
        to_keep = DDE_T == dly_values(iTd);
        tds_matrices{iTd+1} = DDE_B2*diag(to_keep)*DDE_C2;
    end
    
    disp('-- DSA_DDE_chroots: using arnoldi algorithm')
    nb_eig_fctr = DSA.LTI.DDE.nb_eig_fctr;

    % ---------- shifting and scaling ----------
    % a shifted and scaled version of the system is solved. These
    % modifications provide better numerical properties. The resulting
    % characteristic roots are unshifted and rescaled after calculation.
    tds_matrices_ssv = cell(1, nb_matrices); % shifted and scaled version (ssv)
    
    k_shift = -100; % set to 0 to prevent shifting
    max_dly = max(dly_values); % set to 1 to prevent scaling
    
    tds_matrices_ssv{1} = max_dly * (tds_matrices{1} + k_shift * eye(size(tds_matrices{1})));
    for iTd = 1:length(dly_values)
        tds_matrices_ssv{iTd+1} = max_dly * tds_matrices{iTd+1} * exp(k_shift * dly_values(iTd));
    end
    dly_values_ssv = dly_values / max_dly;

    % ---------- system definition ----------
    tds_ssv = tds_create(tds_matrices_ssv, [0, dly_values_ssv]);
    
    % ---------- solving for the roots ----------
    % Initialisation: from the eigenvalues of the no-delay state matrix.
    % Could also be initialised from a random vector.
    nb_extract = ceil(nb_eig_fctr*n*nh_t);
    [all_EIGS_ssv, ~] = tds_arnoldi(tds_ssv, eig(tds_matrices{1}), nb_extract);
    
    % unshifting and unscaling of retrieved eigenvalues
    all_EIGS = all_EIGS_ssv/max_dly - k_shift;

    DSA.LTI.EIGS = all_EIGS;
    
    %% calculation of eigenvectors:
    calc_eigvecs = 1; % calculation of eigenvectors is enabled by default.
    % You may want to deactivate it to accelerate processes, and enable it
    % only when you actually need to check the eigenvectors and to a
    % convergence test.

    convergence_test = DSA.LTI.DDE.convergence_test;
    convergence_eps = DSA.LTI.DDE.convergence_eps;
    if calc_eigvecs
        disp('-- DDE eigs: extraction of eigenvectors')

        nb_extract = length(all_EIGS); % keep this line in case arnoldi == 0
        all_EIGS_converged = zeros(nb_extract, 1);
        all_SVALS = zeros(nb_extract, 1);
        all_SVALS_largest = zeros(nb_extract, 1);
        all_RVEC = zeros(size(DDE_A,1), nb_extract);
        all_LVEC = all_RVEC;

        min_eps = convergence_eps;
        warning('off', 'MATLAB:nearlySingularMatrix')
        fprintf([repmat('.', 1, ceil(nb_extract/10)) '\n'])

        for idx = 1:length(all_EIGS)
            if mod(idx,10) == 1
                fprintf('.')
            end

            lambda = all_EIGS(idx);
            RHS = tds_matrices{1}; % (right hand side)
            for iTd = 1:length(dly_values)
                RHS = RHS + tds_matrices{iTd+1}*exp(-lambda*dly_values(iTd));
            end

            MA = lambda*eye(size(DDE_A)) - RHS;
            [U,S,V] = svds(MA, 1, 'smallest');
            all_SVALS(idx) = S;
            all_LVEC(:, idx) = U;
            all_RVEC(:, idx) = V;
            
            if convergence_test
                [~,S_largest,~] = svds(MA, 1, 'largest','Tolerance',1e-12,'MaxIterations',1e3);
                all_SVALS_largest(idx) = S_largest;
                if S/S_largest <= min_eps
                    all_EIGS_converged(idx) = 1;
                end
            end
        end
        fprintf('\n')
        warning('on', 'MATLAB:nearlySingularMatrix')
        
        DSA.LTI.RVEC = all_RVEC; % columns are the right eigenvectors
        DSA.LTI.LVEC = all_LVEC; % *columns* (not rows) are the left eigenvectors

        if convergence_test
            %%
            % remove the eigenvalues that haven't converged:
            EIGS_conv = DSA.LTI.EIGS;
            EIGS_conv(all_EIGS_converged == 0)    = [];
            
            EIGS_noconv = DSA.LTI.EIGS;
            EIGS_noconv(all_EIGS_converged == 1)    = [];
            
            DSA.LTI.convergence.EIGS_conv = EIGS_conv;
            DSA.LTI.convergence.EIGS_noconv = EIGS_noconv;
            
            % convergence plot
            figure
            grid on
            hold on
            stem(all_SVALS./all_SVALS_largest)
            xlabel('index')
            ylabel('normalised singular value')
            title('Convergence plot')
            
            % plot: identification of converged eigenvalues
            figure
            hold on
            grid on
            plot(real(DSA.LTI.EIGS), imag(DSA.LTI.EIGS)/(2*pi), '.', 'markersize', 12)
            ax = gca;
            ax.ColorOrderIndex = 5;
            plot(real(EIGS_conv), imag(EIGS_conv)/(2*pi), 'o', 'markersize', 4)
            ax.ColorOrderIndex = 2;
            plot(real(EIGS_noconv), imag(EIGS_noconv)/(2*pi), 'o', 'markersize', 4)
            if strcmp(get(groot, 'defaultTextInterpreter'), 'latex')
                xlabel('real part: $\sigma$ [rad/s]', 'fontsize', 14)
                ylabel('imag part: $\omega/(2\pi)$ [Hz]', 'fontsize', 14)
            else
                xlabel('real part: \sigma [rad/s]', 'fontsize', 14)
                ylabel('imag part: \omega/(2\pi) [Hz]', 'fontsize', 14)
            end
            sgtitle('eigenvalues of the system')
        end
    end

    disp('-- DSA_DDE_chroots: done')
end