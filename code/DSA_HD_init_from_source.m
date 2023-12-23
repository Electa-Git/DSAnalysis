function new_HD = DSA_HD_init_from_source(old_HD, new_HD, varnames, new_h_m)
%% new_HD = DSA_HD_init_from_source(old_HD, new_HD, varnames, new_h_m)
% 
% This function writes the harmonic content of variables in varnames from 
% old_HD to new_HD while taking into account the fact that the two  
% structures might have different h_m values.
% 
% Only the variables in varnames are overwritten. The other variables
% in new_HD are left unchanged. old_HD is not changed whatsoever.
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
    
    old_nh_m = length(old_HD.(varnames{1}));
    old_h_m = round((old_nh_m-1)/2);
    old_h2i = @(h) h + old_h_m + 1;
    
    new_nh_m = 2*new_h_m+1;
    new_h2i = @(h) h + new_h_m + 1;
    
    if new_nh_m == old_nh_m
        for vi = 1:length(varnames)
            vname = varnames{vi, 1};
            if isfield(old_HD, vname)
                new_HD.(vname) = old_HD.(vname);
            else
                new_HD.(vname) = 0.01 * ones(new_nh_m,1);
            end
        end
    elseif new_nh_m > old_nh_m
        for vi = 1:length(varnames)
            vname = varnames{vi, 1};
            if isfield(old_HD, vname)
                new_HD.(vname)(new_h2i(-old_h_m:old_h_m)) = old_HD.(vname);
            else
                new_HD.(vname)                            = 0.01 * ones(new_nh_m,1);
            end
            
        end
    else %nh_m_out < nh_m_in
        for vi = 1:length(varnames)
            vname = varnames{vi, 1};
            if isfield(old_HD, vname)
                new_HD.(vname) = old_HD.(vname)(old_h2i(-new_h_m:new_h_m));
            else
                new_HD.(vname) = 0.01 * ones(new_nh_m,1);
            end
        end
    end
    
    disp(['    -> HD_INIT: old h_m = ' num2str(old_h_m) ', new h_m = ' num2str(new_h_m)])
end