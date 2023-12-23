function DSA_save_project(path_inst, varargin)
%% DSA_save_project(path_inst, varargin)
%
% DSA_save_project(path_inst, varargin) is an internal function which you
% normally do not call directly as a user. Use instead:
%   DSA.fcts.save()
%   DSA.fcts.save('example')
%   where the string argument is optional.
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

    if numel(varargin{1}) == 1
        str = ['_' varargin{1}{1}];
    else
        str = '';
    end
    % get a time stamp to avoid overwriting on previous results:
    custom_str = ['_day_' datestr(now,'YY_mm_DD') '_time_' datestr(now,'HH_MM_SS_FFF') str];
    
    % save with detailed file names:
    loc_str = [path_inst '\saved_DSA' custom_str];
    evalin('base',['save(''' loc_str ''', ''DSA'')'])
    
    fprintf('DSA saved at %s\n', loc_str)
end