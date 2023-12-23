function [TA, TB1, TB2, TC1, TC2, TD11, TD12, TD21, TD22, TT_d, TM, TN] = DSA_LTP2DHSS(DSA)
%% [TA, TB1, TB2, TC1, TC2, TD11, TD12, TD21, TD22, TT_d, TM, TN] = DSA_LTP2DHSS(DSA)
%
% This function transforms the generalised LTP system matrices A, B1, B2,
% C1, C2, D11, D12, D21, D22 into their corresponsing toeplitz matrices,
% as defined by the DHSS theory.
% 
% The input matrices are 3D, the third dimension being the harmonic
% dimension. Their depth is equal to 2*h_m+1.
%
% The resulting Toeplitz are defined up to truncation rank h_t.
% The forced periodic rank h_f is also taken into account for the lifting.
% 
% Additionally, time delays provided in vector T_d are taken into account,
% and are included in the lifted representation. Right-multiplication
% delay matrix is provided in the output. Operation must be done as
% follows:
%       TB2  * TM
%       TD12 * TM
%       TD22 * TM
% 
% Time-domain:
% ddt x =  A x +  B1 u +  B2 w
%     y = C1 x + D11 u + D12 w
%     z = C2 x + D21 u + D22 w
%     w = z(t - T)
% 
% Generalised Harmonic State-Space:
% ddt X = (TA-TN) X +  TB1 U +  TB2 * TM W
%     Y =     TC1 X + TD11 U + TD12 * TM W
%     Z =     TC2 X + TD21 U + TD22 * TM W
%     W =     Z(t - TT)
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
d                   = DSA.info.d;
omega_b             = DSA.calc.omega_b;
h_t                 = DSA.calc.h_t;
nh_t                = DSA.calc.nh_t;
h_m                 = DSA.calc.h_m;
h_f                 = DSA.FLIFT.h_f;
LTP_mat             = DSA.LTP.mat;

% -------------------------------------------------------------------------
% From 3D matrices to Toeplitz matrices :

if isempty(h_f)
    % then h_f is automatically set to h_m
    h_f = h_m;
end

TA      = DSA_mat2toeplitz(LTP_mat.A_3D,   h_m, h_f, h_t);
TB1     = DSA_mat2toeplitz(LTP_mat.B1_3D,  h_m, h_f, h_t);
TB2     = DSA_mat2toeplitz(LTP_mat.B2_3D,  h_m, h_f, h_t);
TC1     = DSA_mat2toeplitz(LTP_mat.C1_3D,  h_m, h_f, h_t);
TC2     = DSA_mat2toeplitz(LTP_mat.C2_3D,  h_m, h_f, h_t);
TD11    = DSA_mat2toeplitz(LTP_mat.D11_3D, h_m, h_f, h_t);
TD12    = DSA_mat2toeplitz(LTP_mat.D12_3D, h_m, h_f, h_t);
TD21    = DSA_mat2toeplitz(LTP_mat.D21_3D, h_m, h_f, h_t);
TD22    = DSA_mat2toeplitz(LTP_mat.D22_3D, h_m, h_f, h_t);

% creating matrix N (see LTP2HSS for explanations) :
N_form01 = -h_t:h_t;                        % STEP 01
N_form02 = repmat(N_form01', 1, n)';        % STEP 02
N_form03 = N_form02(:)';                    % STEP 03
TN       = 1j * omega_b * diag(N_form03);   % STEP 04

% building delay matrix TM :
M_form01 = -h_t:h_t;                        % STEP 01
M_form02 = repmat(M_form01', 1, d)';        % STEP 02
M_form03 = M_form02(:)';                    % STEP 03
TT_d = repmat(LTP_mat.T_d, 1, nh_t);
TM = diag(exp(- 1j * omega_b * M_form03 .* TT_d));

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











