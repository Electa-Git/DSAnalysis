function DSA_partfactors(pfa)
%% DSA_partfactors(pfa)
% were pfa stands for participation factors analysis.
% 
% DSA_partfactors calculates participation factors. It is adapted to
% frequency-lifted systems: the names of the variables indicate to which
% time-varying Fourier coefficient they relate.
% 
% Please be aware that not much has been published on participation factors
% of frequency-lifted systems at the time of coding. Make sure you are 
% familiar with the theory of HSS to have some understanding of what 
% frequency-lifted states represent.
% 
% The following fields must be provided in structure pfa.
% EIGS     : eigenvalues
% RVEC     : right eigenvectors
% LVEC     : left eigenvectors
% VARS     : cell of states names
% MODES    : requested modes (first, run with [], then fill in the index of interest)
% XMIN_TXT : indices are displayed for eigenvalues of real part larger than XMIN_TXT
% YMIN_TXT : indices are displayed for eigenvalues of imaginary part larger than YMIN_TXT
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

%% extracting arguments:

EIGS = pfa.EIGS;        % EIGS contains the eigenvalues;
RVEC = pfa.RVEC;        % The columns of RVEC are the right eigenvectors;
LVEC = pfa.LVEC;        % The columns of LVEC are the left  eigenvectors;
% LVEC = (inv(RVEC).'); % is technically another way of getting the same results.
VARS = pfa.VARS;        % VARS contains the names of the states;
MODES = pfa.MODES;      % MODES: empty vector means that all modes are 
                        % displayed. otherwise, give the index of the 
                        % requested modes.
XMIN_TXT = pfa.XMIN_TXT;% Indices are displayed as text on the figure only
                        % for eigenvalues with real part larger than XMIN_TXT.
                        % Setting this to -Inf means that all indices will
                        % be displayed, which can take quite some time.
                        % Set this to a small negative value to display
                        % only the indices of a few modes.
YMIN_TXT = pfa.YMIN_TXT;% Same thing for the y axis, which is convenient 
                        % assuming you have symmetry of modes with respect
                        % to the x axis. Set YMIN_TXT to 0 to have indices
                        % displayed for modes of positive imaginary part
                        % only.

%% Processing MODES

n = length(EIGS);
if isempty(MODES)
    analysed_modes = 1:n;
    detailed_analysis = 1;
else
    analysed_modes = MODES;
    detailed_analysis = 0;
end

cartesian = 1; %#ok<*UNRCH>
% Technically, the code can be extended to the polar case. This is not
% fully implemented yet.

%% figure
if detailed_analysis
    figure
    hold on
    grid on
    if cartesian
        plot(real(EIGS), imag(EIGS)/(2*pi), '.', 'markersize', 12) % normalisation by 2*pi to have Hz on the imaginary axis
    else
        plot(real(EIGS), imag(EIGS), '.', 'markersize', 12)
    end
    for idx = 1:length(EIGS)
        ev = EIGS(idx);
        if real(ev) > XMIN_TXT && imag(ev)/(2*pi) > YMIN_TXT
            text(double(real(ev)), double(imag(ev)/(2*pi)), sprintf('%04d',idx), 'fontsize', 11)
        end
    end
    set(gca,'ColorOrderIndex', 2)
    if cartesian
        plot([0, 0], ylim, '--', 'LineWidth', 1.5);
    else
        th = 0:pi/100:2*pi; xunit = cos(th); yunit = sin(th); plot(xunit, yunit, '--', 'LineWidth', 1.5);
    end

    if strcmp(get(groot, 'defaultTextInterpreter'), 'latex')
        xlabel('real part: $\sigma$ [rad/s]', 'fontsize', 14)
        ylabel('imag part: $\omega/(2\pi)$ [Hz]', 'fontsize', 14)
    else
        xlabel('real part: \sigma [rad/s]', 'fontsize', 14)
        ylabel('imag part: \omega/(2\pi) [Hz]', 'fontsize', 14)
    end
    sgtitle('eigenvalues/characteristic roots')
end

%% sorting modes
if detailed_analysis
    sorting_matrix = zeros(n, 2);
    sorting_matrix(:, 1) = EIGS;
    sorting_matrix(:, 2) = 1:n;

    if cartesian % sorting based on real and imaginary parts
        sorting_matrix(:, 3) = real(EIGS);
        sorting_matrix(:, 4) = imag(EIGS);
        [sorted_matrix, ~] = sortrows(sorting_matrix, [3 4]);
        sorted_eig = sorted_matrix(:, 1);
        sorted_idx = sorted_matrix(:, 2);

        for idx = 1:n
            ev = sorted_eig(idx);
            ei = sorted_idx(idx);

            if isreal(ev)
                fprintf('mode %5d : real part: %12.4f (real-valued)\n', ei, ev)
            else
                ev_freq = imag(ev) / (2*pi);
                ev_damp = - real(ev) / sqrt(real(ev)^2 + imag(ev)^2);
                fprintf('mode %5d : real part: %12.4f,    freq: %10.2f Hz,    damping: %10.5f\n', ei, real(ev), ev_freq, ev_damp)
            end
        end
    else % sorting based on magnitude and argument
        sorting_matrix(:, 3) = abs(EIGS);
        sorting_matrix(:, 4) = arg(EIGS);
        [sorted_matrix, ~] = sortrows(sorting_matrix, [3 4]);
        sorted_idx = sorted_matrix(:, 2);
        sorted_mag = sorted_matrix(:, 3);
        sorted_arg = sorted_matrix(:, 4);

        for idx = 1:n
            ei = sorted_idx(idx);
            em = sorted_mag(idx);
            ea = sorted_arg(idx);

            fprintf('eigenvalue %5d : mag: %12.4f,    arg: %10.2f\n', ei, em, ea)
        end
    end
end
% ---------------------------------------------------------------------
% participation analysis

% The left eigenvectors carry information about the CONTROLLABILITY of
% individual modal variable by individual state variables. If the
% eigenvectors are normalized then the kth element L_ki of the ith left
% eigenvector determines the magnitude and phase of the share of the kth 
% state variable x_k(t) in the activity of the ith mode z_i(t).

% The right eigenvectors carry information about the OBSERVABILITY of 
% individual modal variables in individual state variables. If the 
% eigenvectors are normalized then R_ki determines the magnitude and phase 
% of the share of the ith modal variable z_i(t) in the activity of the kth 
% state variable x_k(t)

% Coefficients p_ki = L_ki*R_ki are referred to as the participation
% factors. Each participation factor is a product of the kth element of the
% ith left and right eigenvectors. It quantifies the sensitivity of the ith
% eigenvalue to the kth diagonal element of the state matrix. Element R_ki
% contains information about the observability of the ith modal variable in
% the kth state variable, while L_ki contains information about the
% controllability of the ith modal variable using the kth state variable.
% Hence the product p_ki = L_ki*R_ki contains information about the
% observability and controllability. Consequently, the participation factor 
% is a good measure of correlation between the ith modal variable and the kth 
% state variable.

P = LVEC.*RVEC;
P_abs = abs(P);
% In P, columns correspond to modes while rows correspond to states.
% In a given row/for a given state, a larger value indicates which mode/col
% is more correlated to the state.

% The right and left eigenvectors are already orthogonal so scalar product
% of LVEC(i) and RVEC(j) is zero if i~=j. However, we also want this scalar
% product to be equal to 1 when i=j. To do so, each vector is normalised by
% the sqrt of the value of the scalar product. Equivalently, each vector of
% the matrix of participation factors can be normalised by the value of the
% scalar product directly.

% Normalisation for each column (i.e. each mode):
P_norm = zeros(size(P_abs));
for idx = 1:n
    P_norm(:, idx) = P_abs(:, idx) / sum(P_abs(:, idx));
end

%%
% ---------------------------------------------------------------------
% Tables of analysis
max_cumulative_sum = 0.99; % parameter for the display of the Tables.

fprintf('\n\nParticipation factors analysis based on magnitude of elements.\n')
fprintf('X%% = correlation between the state and the mode.\n')
fprintf('  mode | states --\n')
fprintf('   |\n')
for mode = analysed_modes
    [sort_B, sort_I] = sort(P_norm(:, mode), 'descend');
    cumulative_sum_sort_B = cumsum(sort_B);
    until = find(cumulative_sum_sort_B > max_cumulative_sum, 1);

    if ~isempty(until)
        maxlength = max(cellfun('length', VARS(sort_I(1:until)))); % max length of displayed names
        dispval = max(maxlength, 15); % increase here if the length of the variables names is generally larger than the default value
        
        fprintf([' %3d     ' repmat([' %' num2str(dispval) 'd'      ], 1, until) '\n'], mode, sort_I(1:until))              % the states indices
        fprintf(['         ' repmat([' %' num2str(dispval) 's'      ], 1, until) '\n'], VARS{sort_I(1:until)})              % the states names
        fprintf([' correl  ' repmat([' %' num2str(dispval-1) '.1f%%'], 1, until) '\n'], 100*sort_B(1:until))                % the participation
        fprintf([' cumul   ' repmat([' %' num2str(dispval-1) '.1f%%'], 1, until) '\n'], 100*cumulative_sum_sort_B(1:until)) % the cumulative sum
        fprintf(['_________' repmat(repmat('_', 1, dispval+1), 1, until) '\n'])
    else % blanks(4)
        fprintf(' %3d           all \n', mode)
        fprintf(' corr          0.0%% \n')
        fprintf(' cumul          / \n')
        fprintf('__________________\n')
    end
end
end