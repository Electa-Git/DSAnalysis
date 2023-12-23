function DSA = DSA_sequences(DSA, HD_source, TD_source, vnames)
%% DSA = DSA_sequences(DSA, HD_source, TD_source, vnames)
%
% This function relies on Fortescue transformation to calculate the
% positive-, negative- and zero-sequence components of three-phase "abc"
% sets of variables. The harmonic content of the three sequences is
% calculated first in the harmonic domain, relying on the HD data in
% HD_source. The results are stored in the same structure.
%
% Next, the corresponding TD waveforms are obtained via evaluation of the
% Fourier-series formula via DSA_HD2TD(). The waveforms are stored in
% TD_source.
%
% Both HD_source and TD_source are string names which must exist in DSA.HD
% and DSA.TD beforehand.
% Cell vnames lists the name of the variables to be processed. For example,
% vnames = {'v_g';
%           'i_s'};
% calculates the sequences of v_g_a, v_g_b, v_g_c, and then the sequences
% of i_s_a, i_s_b, i_s_c. The code currently only works for variables of
% which the _a, _b and _c subscript appear at the end of the variable
% names. For this example, the sequences are stored in v_g_ps, v_g_ns,
% v_g_zs, and then in i_s_ps, i_s_ns, i_s_zs.
%
% Extension:
% Use vnames = {{'v_g', 'ref'};
%               'i_s'};
% to make a distinction between what comes before and after the phase 
% letter. For this example, the code would calculate the sequences of
% v_g_a_ref, v_g_b_ref, v_g_c_ref and then those of 
% i_s_a, i_s_b, i_s_c.
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

disp('--- sequence calculations:')

sigma = exp(1j*2*pi/3);
% DIFFERENT FORTESCUE TRANSFORMATIONS FOR POSITIVE AND NEGATIVE FREQUENCIES!
FORTESCUE_pos_freq = 1/3 * [1, 1, 1;
           1, sigma  , sigma^2;
           1, sigma^2, sigma];
FORTESCUE_neg_freq = 1/3 * [1, 1, 1;
           1, sigma^2, sigma;
           1, sigma  , sigma^2];

HD_src = DSA.HD.(HD_source);
TD_src = DSA.TD.(TD_source);

for idx = 1:size(vnames, 1)
    vname_line = vnames{idx,:};
    if ischar(vname_line)
        vname = vname_line;
        postfix = '';
    else
        vname = vname_line{1};
        postfix = ['_' vname_line{2}];
    end

    vec_abc = [(HD_src.([vname '_a' postfix])).';
               (HD_src.([vname '_b' postfix])).';
               (HD_src.([vname '_c' postfix])).'];
    
    vec_pnz = [FORTESCUE_neg_freq * vec_abc(:,(-HD_src.h_m:-1) + HD_src.h_m+1),...
               FORTESCUE_pos_freq * vec_abc(:, (0:HD_src.h_m)  + HD_src.h_m+1)];
    
    HD_src.([vname '_zs' postfix]) = vec_pnz(1,:).';
    HD_src.([vname '_ps' postfix]) = vec_pnz(2,:).';
    HD_src.([vname '_ns' postfix]) = vec_pnz(3,:).';
    
    TD_src.([vname '_zs' postfix]) = DSA_HD2TD(HD_src.([vname '_zs' postfix]), TD_src.t, HD_src.f_b);
    TD_src.([vname '_ps' postfix]) = DSA_HD2TD(HD_src.([vname '_ps' postfix]), TD_src.t, HD_src.f_b);
    TD_src.([vname '_ns' postfix]) = DSA_HD2TD(HD_src.([vname '_ns' postfix]), TD_src.t, HD_src.f_b);
    
	if ~any(strcmp(mycell, [vname '_zs' postfix])) % checks if the variables are already in the cell
		DSA.vars.lasi_outputs = [DSA.vars.lasi_outputs;
								[vname '_zs' postfix];
								[vname '_ps' postfix];
								[vname '_ns' postfix]];
	end
    
    disp([vname '_a' postfix ']    [' vname '_zs' postfix])
    disp([vname '_b' postfix '| -> |' vname '_ps' postfix])
    disp([vname '_c' postfix ']    [' vname '_ns' postfix])
    disp('')
end

DSA.HD.(HD_source) = HD_src;
DSA.TD.(TD_source) = TD_src;

disp('--- sequence calculations: done')
end