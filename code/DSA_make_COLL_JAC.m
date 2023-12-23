function DSA_make_COLL_JAC(DSA)
%% DSA_make_COLL_JAC(DSA)
%
% This function makes the analytical Jacobian of the Fourier-based
% collocation formulation.
% Specifically, it relies on the symbolic nonlinear system to obtain the
% derivatives with respect to the sates and to-delay variables, if any.
%
% Author: Philippe De Rua
% Updated: 2023-06-06

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

sstm = DSA.sstm;
inst = DSA.inst;

if DSA.data.delay.is_delayed
    function_name_short = sprintf('%s_%s_COLL_JAC_DLYfun',sstm,inst);
    function_name       = sprintf('sstms/%s/%s/%s_%s_COLL_JAC_DLYfun.m',sstm,inst,sstm,inst);
    script_name         = sprintf('sstms/%s/%s/%s_%s_COLL_JAC_DLY.m',sstm,inst,sstm,inst);
else
    function_name_short = sprintf('%s_%s_COLL_JACfun',sstm,inst);
    function_name       = sprintf('sstms/%s/%s/%s_%s_COLL_JACfun.m',sstm,inst,sstm,inst);
    script_name         = sprintf('sstms/%s/%s/%s_%s_COLL_JAC.m',sstm,inst,sstm,inst);
end

% retrieve expression of f(x,u), such that the COLL solves dxdt = f(x,u).
% if there are delayed variables, the vector is extended with delayed
% variables:
if DSA.data.delay.is_delayed
    COLL_dxdt = [DSA.NLSS.dxdt;
                 DSA.NLSS.lasi.z];
    COLL_unknowns = [DSA.vars.lasi_states;
                     DSA.vars.lasi_delayed];
else
    COLL_dxdt = DSA.NLSS.dxdt;
    COLL_unknowns = DSA.vars.lasi_states;
end

% get symbolic jacobian of the COLL problem, i.e. dF(x)/dx
disp('--- COLL_JAC: symbolic calculation of COLL jacobian')
J_symb = jacobian(COLL_dxdt, cell2sym(COLL_unknowns));
disp('--- COLL_JAC: symbolic COLL jacobian calculated')

% get a "zero symbolic object":
% this is necessary to pad the matlab function and dynamically extend to
% numerical matrix with arbitrary maximum harmonic ranks when solving the
% collocation system
syms null_symb;
Z_symb = null_symb*ones(size(J_symb));

% EXAMPLE
% Let us consider a simple case without delay: the linearised system is:
% A(t)x(t) =
% [A11(t) A12(t)] [x1(t)]
% [A21(t) A22(t)] [x2(t)]
% The collocation method is such that each entry is evaluated over a grid
% of points [t1, t2, t3].
% In this code, the expansion is dynamic. There are basically two main
% "padding" options, but from an implementation-perspective, the easiest is
% to "keep variables by themselves", instead of keeping "time points" by
% themselves:
% The sampled-time system becomes:
%
% [A11(t1)              A12(t1)              ] [x1(t1)]
% [       A11(t2)              A12(t2)       ] [x1(t2)]
% [              A11(t3)              A12(t3)] [x1(t3)]
% [A21(t1)              A22(t1)              ] [x2(t1)]
% [       A21(t2)              A22(t2)       ] [x2(t2)]
% [              A21(t3)              A22(t3)] [x2(t3)]

% Convert the symbolic jacobian into a matlab function (rewritten in a
% script for simplicity)
% this may take some time. In theory, this matlab function is created once  
% every time the differential equations are modified.

disp('--- COLL_JAC: creation of COLL JAC function (this may take some time)')
matlabFunction(J_symb + Z_symb, 'file', function_name, 'optimize', false);
disp('--- COLL_JAC: COLL JAC function created')

% editing the matlabfunction ----------------------------------------------
% --- Open file with fopen()
fileID_in = fopen(function_name, 'rt');
fileID_out = fopen(script_name, 'w');

% --- read first line of file
textLine = fgetl(fileID_in);
functionCallLine = textLine;

out = regexp(functionCallLine, '[^()]*','match');
fprintf(fileID_out, 'J_RHS = %s(%s);\n\n', function_name_short, out{2});

% --- run
line_nb = 1;
time_lifting_done = 0;

while ischar(textLine) && ~startsWith(strtrim(textLine), 'out1 = reshape')
     % --> COPY EVERYTHING UNTIL THE VERY LAST OUT1 IS FOUND

    if startsWith(strtrim(textLine), 'mt') && ~contains(textLine, '[mt')
        % in this case, the textline must be modified, otherwise not.
        old_str = {'['; ','; ']'};
        new_str = {'[{diag('; ')},{diag('; ')}]'};
        textLine = replace(textLine, old_str, new_str);
        time_lifting_done = 1;
    end

    fprintf(fileID_out, '%s\n', textLine);
    textLine = fgetl(fileID_in);
    line_nb = line_nb + 1;
end

% --- Dealing with the reshape function
% (making a cellarray, reshaping the cellarray, transforming the cellarray into a matrix)

% recognise the main array and the size array
out = regexp(textLine, '[^[]]*','match');
cellarray_str = ['[' out{2} ']'];
out = regexp(textLine, '[^(,);]*','match');
sz1 = out{length(out)-1};
sz2 = out{length(out)};

if ~time_lifting_done
    % for very simple systems, there may not be intermediate mt matrices.
    % now is thus the time to time-lift to sample points:
    old_str = {'['; ','; ']'};
    new_str = {'[{diag('; ')},{diag('; ')}]'};
    cellarray_str = replace(cellarray_str, old_str, new_str);
end
fprintf(fileID_out, 'cellarray = %s;\n', cellarray_str);
fprintf(fileID_out, 'out1 = cell2mat(reshape(cellarray,%s,%s));\n', sz1, sz2);
fprintf(fileID_out, 'end\n');

% --- Close both files with fclose()
fclose(fileID_in);
fclose(fileID_out);
% --- Delete the input file with delete().
delete(function_name);

disp('--- COLL_JAC: analytical collocation jacobian written at:')
disp(['            ' script_name])
% ---------------------------------------------------------- end of edition
end