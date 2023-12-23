function TM = DSA_mat2toeplitz(M, h_m, h_f, h_t)
%% TM = DSA_mat2toeplitz(M, h_m, h_f, h_t)
%
% DSA_mat2toeplitz is a function at the core of the frequency-lifting
% method: it provides the block-Toeplitz form of a matrix based on its set
% of matrix Fourier coefficients.
%
% M is a 3D matrix, and it corresponds to any of A, B, C or D in a
% state-space representation.
% The first and second dimensions of the matrix are related to the system
% dimensions, i.e. numbers of inputs and/or outputs and/or states.
% The 3rd dimension of 3D matrix M is the harmonic dimension:
% 
% The 3D matrices can be seen as stacks of harmonic layers. There are
% 2*h_m+1 layers. The first operation is to take the layers and make a 
% series out of them. The series of layers is stored in a cell. Each 
% element of the cell named M_series is therefore a matrix.
% 
% Due to the mapping between series and matrix configurations, M_series 
% actually has 4*h_t+1 elements.
% 
% IMPORTANT:
% If h_t >= h_m, then the 4*h_t+1 elements are obtained from the 2*h_m+1
% elements of M, and extra spaces are filled with zero matrices.
% If h_t < h_m, there are two possibilities:
%       - strict_truncation == 1 will make that extra elements in M will be
%         ignored, and replaced by zero matrices
%       - strict_truncation == 0 will keep the extra elements in M, and
%         also fill remaining space with zero matrices
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

% basicaly set h_f automatically to either h_m or h_t
if h_f < h_m
    M = M(:,:, h_m+1+(-h_f:h_f));
    h_m = h_f;
end

% -------------------------------------------------------------------------
% Easiest case when the harmonic matrix only contains one harmonic layer:
% It is not a 3D matrix but a 2D matrix.

if h_m == 0
    TM = kron(eye(2*h_t+1), M);

% -------------------------------------------------------------------------
% Otherwise it is a 3D matrix:
else
    % (1) initialise series-cell with zero matrices
    
    [nb_rows, nb_cols] = size(M(:,:,1));
    M_series = repmat({zeros(nb_rows, nb_cols)}, 4*h_t+1, 1); 
    
    % (2) determine which harmonics of M must be copied into M_series
    
    if h_t >= h_m
        h_limit = h_m;
    else
        h_limit = min(h_m, 2*h_t);
    end
    
    % (3) copy from M to M_series accounting for the harmonic index
    % mappings:
    
    for h = -h_limit:h_limit
        M_series{2*h_t+1+h} = M(:,:,h_m+1+h);
    end
    
    % (4) Apply mapping between M_series and TM_cell, then transform 
    % TM_cell into a matrix:
    %
    % Example of mapping for h_m = h_t = 2:
    % 
    % harmonic_map = 
    %      [0   -1   -2    0    0
    %       1    0	 -1   -2    0
    %       2    1	  0   -1   -2
    %       0    2	  1    0    1
    %       0    0	  2    1    0]
    % 
    % index_map =
    %      [5    4    3    2	1
    %       6    5    4    3    2
    %       7    6    5    4    3
    %       8    7    6    5    4
    %       9    8    7    6    5]
    
    index_map = toeplitz((2*h_t+1):(4*h_t+1), (2*h_t+1):-1:1);
    TM_cell = M_series(index_map);
    TM = cell2mat(TM_cell);
end
end