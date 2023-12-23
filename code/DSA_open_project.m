function DSA_open_project(str_DATA, str_VARS, str_HD, str_uOft, str_EQNS, varargin)
%% DSA_open_project(str_DATA, str_VARS, str_HD, str_uOft, str_EQNS, varargin)
%
% DSA_open_project is an internal function which you
% normally do not call directly as a user. Use instead:
%   DSA.fcts.open()
%   DSA.fcts.open('DATA')
%   DSA.fcts.open('EQNS')
%   DSA.fcts.open('VARS')
%   DSA.fcts.open('HD')
%   DSA.fcts.open('uOft')
%   where the string arguments are optional. If no argument is provided,
%   all files related to the project described in DSA will be opened.
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

    if numel(varargin{:}) ~= 1
        if isfile(str_DATA)
             open(str_DATA);
        end
        if isfile(str_VARS)
             open(str_VARS);
        end
        if isfile(str_HD)
             open(str_HD);
        end
        if isfile(str_uOft)
             open(str_uOft);
        end
        if isfile(str_EQNS)
             open(str_EQNS);
        end
    else
        if     strcmp(varargin{1}{1}, 'DATA') && isfile(str_DATA)
             open(str_DATA);
        elseif strcmp(varargin{1}{1}, 'VARS') && isfile(str_VARS)
             open(str_VARS);
        elseif strcmp(varargin{1}{1}, 'HD')   && isfile(str_HD)
             open(str_HD);
        elseif strcmp(varargin{1}{1}, 'uOft') && isfile(str_uOft)
             open(str_uOft);
        elseif strcmp(varargin{1}{1}, 'EQNS') && isfile(str_EQNS)
             open(str_EQNS);
        end
    end
end