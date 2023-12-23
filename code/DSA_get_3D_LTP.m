function DSA = DSA_get_3D_LTP(DSA)
%% DSA = DSA_get_3D_LTP(DSA)
%
% This function performs a Fourier analysis of the time-periodic matrix
% coefficients of the linearised state-space system.
% 1) All variables are evaluated as numerical sampled-time vectors
% 2) All parameters are evaluated as numerical values
% 3) The symbolic state-space coefficients are evaluated numerically and
% element-wise for each matrix: each matrix entry becomes a sampled-time
% vector.
% 4) The DFT (discrete fourier transform) is applied to each sampled-time
% vector in order to obtain its harmonic spectrum.
% 5) The harmonic spectra are stored in a 3D matrix: the first and second
% dimensions are those of the original matrix coefficients. The third
% dimension is the harmonic index dimension, from - to + frequencies.
% 6) The resulting matrix contains information of a linear time-periodic
% (LTP) system, but stored as constant data in a 3D matrix, hence the name
% of the function.
% 
% Evaluating every matrix entry element-wise is highly convenient and
% straightforward approach. However, it is not very efficient: this
% function is not optimised for speed.
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

%#ok<*NASGU>

DSA = DSA_set(DSA); % updating DSA in case numerical parameter values have been changed.

nh_m    = DSA.calc.nh_m;
TD2HD   = DSA.calc.TD2HD;
HD2TD   = DSA.calc.HD2TD;
n       = DSA.info.n;
m       = DSA.info.m;
p       = DSA.info.p;
d       = DSA.info.d;
LNSS    = DSA.LNSS;
h2i_h_m = DSA.fcts.h2i_h_m;

%% preparing variables and parameters for evaluation

% ---------------------------------------------------------------------
% variables are assigned with their sampled-time vector:

if isfield(DSA.HD, DSA.LTP.HD_source)
    % retrieving HD source:
    HD_struct = DSA.HD.(DSA.LTP.HD_source);

    % evaluating variables numerically as sampled-time vectors:
    if DSA.data.delay.is_delayed
        varnames = [DSA.vars.lasi_inputs;
                    DSA.vars.lasi_states;
                    DSA.vars.lasi_delayed];
    else
        varnames = [DSA.vars.lasi_inputs;
                    DSA.vars.lasi_states];
    end

    for idx = 1:length(varnames)
        vname = varnames{idx};
        eval([vname ' = real(HD2TD*HD_struct.(vname));']);
    end

    disp(['--- DSA_get_3D_LTP: evaluation along the periodic trajectory defined in HD_source: ' DSA.LTP.HD_source])
else
    warning('DSA_get_3D_LTP: Requested HD_source is not available in the data structure. This is fine for LTI systems, but will most probably lead to an error for LTP and NTIp systems.')
end

% ---------------------------------------------------------------------
% parameters are assigned with their numerical value:

varnames = fields(DSA.data.param);
for idx = 1:length(varnames)
    vname = varnames{idx};
    eval([vname ' = DSA.data.param.(vname);']);
end

%% From symbolic to harmonic matrices

disp('--- DSA_get_3D_LTP: from 2D TD symbolic to 3D HD matrices')

% ---------------------------------------------------------------------
% initializing empty 3D matrices, distinguishing between delayed and
% non-delayed cases.

if ~DSA.data.delay.is_delayed
    A_3D = zeros(n, n, nh_m);
    B_3D = zeros(n, m, nh_m);
    C_3D = zeros(p, n, nh_m);
    D_3D = zeros(p, m, nh_m);
else
    A_3D   = zeros(n, n, nh_m);
    B1_3D  = zeros(n, m, nh_m);
    B2_3D  = zeros(n, d, nh_m);
    C1_3D  = zeros(p, n, nh_m);
    C2_3D  = zeros(d, n, nh_m);
    D11_3D = zeros(p, m, nh_m);
    D12_3D = zeros(p, d, nh_m);
    D21_3D = zeros(d, m, nh_m);
    D22_3D = zeros(d, d, nh_m);
    T_d      = DSA.data.delay.T_d;
end

% ---------------------------------------------------------------------
% filling in matrices with the harmonic data

% Each expression of the symbolic matrices is evaluated in TD with
% sample-points. The symbolic expressions are vectorized and then evaluated
% numerically.

if ~DSA.data.delay.is_delayed
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - A_3D
    fprintf('    calculating A(t)')
    for r = 1:n
        for c = 1:n
            TD_val = eval(vectorize(LNSS.mat.A(r, c)));
            if isscalar(TD_val) % constant value, no need for FFT
                A_3D(r, c, h2i_h_m(0)) = TD_val;
            else
                A_3D(r, c, :) = TD2HD * TD_val;
            end
        end
        if mod(r,10)==0
            fprintf('.')
        end
    end
    fprintf('\n')
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - B_3D
    disp('    calculating B(t)')
    for r = 1:n
        for c = 1:m
            TD_val = eval(vectorize(LNSS.mat.B(r, c)));
            if isscalar(TD_val)
                B_3D(r, c, h2i_h_m(0)) = TD_val;
            else
                B_3D(r, c, :) = TD2HD * TD_val;
            end
        end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - C_3D
    disp('    calculating C(t)')
    for r = 1:p
        for c = 1:n
            TD_val = eval(vectorize(LNSS.mat.C(r, c)));
            if isscalar(TD_val)
                C_3D(r, c, h2i_h_m(0)) = TD_val;
            else
                C_3D(r, c, :) = TD2HD * TD_val;
            end
        end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - D_3D
    disp('    calculating D(t)')
    for r = 1:p
        for c = 1:m
            TD_val = eval(vectorize(LNSS.mat.D(r, c)));
            if isscalar(TD_val)
                D_3D(r, c, h2i_h_m(0)) = TD_val;
            else
                D_3D(r, c, :) = TD2HD * TD_val;
            end
        end
    end
    
    DSA.LTP.mat.A_3D = A_3D;
    DSA.LTP.mat.B_3D = B_3D;
    DSA.LTP.mat.C_3D = C_3D;
    DSA.LTP.mat.D_3D = D_3D;

else
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - A_3D
    fprintf('   calculating A(t)')
    for r = 1:n
        for c = 1:n
            TD_val = eval(vectorize(LNSS.mat.A(r, c)));
            if isscalar(TD_val)
                A_3D(r, c, h2i_h_m(0)) = TD_val;
            else
                A_3D(r, c, :) = TD2HD * TD_val;
            end
        end
        if mod(r,10)==0
            fprintf('.')
        end
    end
    fprintf('\n')
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - B1_3D
    disp('   calculating B1(t)')
    for r = 1:n
        for c = 1:m
            TD_val = eval(vectorize(LNSS.mat.B1(r, c)));
            if isscalar(TD_val)
                B1_3D(r, c, h2i_h_m(0)) = TD_val;
            else
                B1_3D(r, c, :) = TD2HD * TD_val;
            end
        end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - B2_3D
    disp('   calculating B2(t)')
    for r = 1:n
        for c = 1:d
            TD_val = eval(vectorize(LNSS.mat.B2(r, c)));
            if isscalar(TD_val)
                B2_3D(r, c, h2i_h_m(0)) = TD_val;
            else
                B2_3D(r, c, :) = TD2HD * TD_val;
            end
        end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - C1_3D
    disp('   calculating C1(t)')
    for r = 1:p
        for c = 1:n
            TD_val = eval(vectorize(LNSS.mat.C1(r, c)));
            if isscalar(TD_val)
                C1_3D(r, c, h2i_h_m(0)) = TD_val;
            else
                C1_3D(r, c, :) = TD2HD * TD_val;
            end
        end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - C2_3D
    disp('   calculating C2(t)')
    for r = 1:d
        for c = 1:n
            TD_val = eval(vectorize(LNSS.mat.C2(r, c)));
            if isscalar(TD_val)
                C2_3D(r, c, h2i_h_m(0)) = TD_val;
            else
                C2_3D(r, c, :) = TD2HD * TD_val;
            end
        end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - D11_3D
    disp('   calculating D11(t)')
    for r = 1:p
        for c = 1:m
            TD_val = eval(vectorize(LNSS.mat.D11(r, c)));
            if isscalar(TD_val)
                D11_3D(r, c, h2i_h_m(0)) = TD_val;
            else
                D11_3D(r, c, :) = TD2HD * TD_val;
            end
        end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - D12_3D
    disp('   calculating D12(t)')
    for r = 1:p
        for c = 1:d
            TD_val = eval(vectorize(LNSS.mat.D12(r, c)));
            if isscalar(TD_val)
                D12_3D(r, c, h2i_h_m(0)) = TD_val;
            else
                D12_3D(r, c, :) = TD2HD * TD_val;
            end
        end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - D21_3D
    disp('   calculating D21(t)')
    for r = 1:d
        for c = 1:m
            TD_val = eval(vectorize(LNSS.mat.D21(r, c)));
            if isscalar(TD_val)
                D21_3D(r, c, h2i_h_m(0)) = TD_val;
            else
                D21_3D(r, c, :) = TD2HD * TD_val;
            end
        end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - D22_3D
    disp('   calculating D22(t)')
    for r = 1:d
        for c = 1:d
            TD_val = eval(vectorize(LNSS.mat.D22(r, c)));
            if isscalar(TD_val)
                D22_3D(r, c, h2i_h_m(0)) = TD_val;
            else
                D22_3D(r, c, :) = TD2HD * TD_val;
            end
        end
    end

    DSA.LTP.mat.A_3D    = A_3D;
    DSA.LTP.mat.B1_3D   = B1_3D;
    DSA.LTP.mat.B2_3D   = B2_3D;
    DSA.LTP.mat.C1_3D   = C1_3D;
    DSA.LTP.mat.C2_3D   = C2_3D;
    DSA.LTP.mat.D11_3D  = D11_3D;
    DSA.LTP.mat.D12_3D  = D12_3D;
    DSA.LTP.mat.D21_3D  = D21_3D;
    DSA.LTP.mat.D22_3D  = D22_3D;
    DSA.LTP.mat.T_d     = T_d;
end

disp('--- DSA_get_3D_LTP: 3D harmonic matrices defined')

%% Spectral norm (maximum singular value) of matrix Fourier coefficients:
% the calculation is done only for the state matrix A(t), however it is
% easily extended to other state-space matrices in both delayed and
% non-delayed representations.

A_3D_normalised = A_3D / svds(A_3D(:,:,h2i_h_m(0)),1);
fprintf('\nNorm (maximum singular value) of matrix Fourier coefficients of A(t):\n')
fprintf('    index, norm\n');
for idx = 0:DSA.calc.h_m
    val = svds(A_3D_normalised(:,:,DSA.fcts.h2i_h_m(idx)),1);
    if val < 1e-15
        fprintf('    %2d,    %d\n', idx, 0);
    else
        fprintf('    %2d,    %.2e\n', idx, val);
    end
end

disp('--- DSA_get_3D_LTP: done')
end













