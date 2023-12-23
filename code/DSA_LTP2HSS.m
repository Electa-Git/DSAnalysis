 function [TA, TB, TC, TD, TN] = DSA_LTP2HSS(DSA)
%% [TA, TB, TC, TD, TN] = DSA_LTP2HSS(DSA)
%
% This function transforms the LTP system matrices A, B, C, D into their
% corresponsing Toeplitz matrices, as defined by the HSS theory.
%
% The input matrices are 3D, the third dimension being the harmonic
% dimension. Their depth is equal to 2*h_m+1.
% 
% The resulting Toeplitz are defined up to truncation rank h_t.
% The forced periodic rank h_f is also taken into account for the lifting.
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

n                   = DSA.info.n;
omega_b             = DSA.calc.omega_b;
h_t                 = DSA.calc.h_t;
h_m                 = DSA.calc.h_m;
h_f                 = DSA.FLIFT.h_f;
LTP_mat             = DSA.LTP.mat;

% -------------------------------------------------------------------------
% From 3D matrices to Toeplitz matrices :

if isempty(h_f)
    % then h_f is automatically set to h_m
    h_f = h_m;
end

TA = DSA_mat2toeplitz(LTP_mat.A_3D, h_m, h_f, h_t);
TB = DSA_mat2toeplitz(LTP_mat.B_3D, h_m, h_f, h_t);
TC = DSA_mat2toeplitz(LTP_mat.C_3D, h_m, h_f, h_t);
TD = DSA_mat2toeplitz(LTP_mat.D_3D, h_m, h_f, h_t);

N_form01 = -h_t:h_t;                    % STEP 01 (see explanation below)
N_form02 = repmat(N_form01', 1, n)';    % STEP 02
N_form03 = N_form02(:)';                % STEP 03
TN = 1j * omega_b * diag(N_form03);     % STEP 04

% -------------------------------------------------------------------------
% How matrix TN is obtained: with a series of MatLab tricks:
% (example with h_t = 2 and n = 3):
% 
% STEP 01:
% --------
%    [-2    -1     0     1     2]
% 
% STEP 02:
% --------
%    [-2    -1     0     1     2
%     -2    -1     0     1     2
%     -2    -1     0     1     2]
% 
% STEP 03:
% --------
%  [-2    -2    -2    -1    -1    -1     0     0     0      1     1     1     2     2     2]
% 
% STEP 04 (without 1j*omega_1):
% -------------------------------
%    [-2     0     0     0     0     0     0     0     0     0     0     0     0     0     0
%      0    -2     0     0     0     0     0     0     0     0     0     0     0     0     0
%      0     0    -2     0     0     0     0     0     0     0     0     0     0     0     0
%      0     0     0    -1     0     0     0     0     0     0     0     0     0     0     0
%      0     0     0     0    -1     0     0     0     0     0     0     0     0     0     0
%      0     0     0     0     0    -1     0     0     0     0     0     0     0     0     0
%      0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
%      0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
%      0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
%      0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
%      0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
%      0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
%      0     0     0     0     0     0     0     0     0     0     0     0     2     0     0
%      0     0     0     0     0     0     0     0     0     0     0     0     0     2     0
%      0     0     0     0     0     0     0     0     0     0     0     0     0     0     2]
end