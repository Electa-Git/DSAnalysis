function status = DSA_odeprog(t, y, flag, varargin) %#ok<INUSL>
%% status = DSA_odeprog(t, y, flag, varargin)
%
% Function DSA_odeprog helps keeping track of the progress of ODE solvers.
% 
% Original authors:
% Tim Franklin
% Virginia Tech
% Jesse Norris
% Wake Forrest
% May 2006
% 
% --> updated by Philippe De Rua in 2022
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

global odeprogglobvar

if nargin < 3 || isempty(flag) 
    if etime(clock,odeprogglobvar(9:14)) > 2 % update every two seconds
        t_start = odeprogglobvar(1);
        t_stop  = odeprogglobvar(2);
        sstrt   = odeprogglobvar(3:8);
        perc    = (t(end) - t_start)/(t_stop-t_start);
        
        fprintf('%08.4f%% -- %s -- %s\n', perc*100, etimev(clock,sstrt), etimev(etime(clock,sstrt)/perc*(1-perc)));
        
        odeprogglobvar(9:14)=clock;
    end
else
    switch(flag)
    case 'init'  
        odeprogglobvar       = zeros(1,14);
        odeprogglobvar(1)    = t(1);
        odeprogglobvar(2)    = t(end);
        odeprogglobvar(3:8)  = clock;
        odeprogglobvar(9:14) = clock;
        sstrt = odeprogglobvar(3:8);
        perc = 0;
        fprintf('%s -- %s -- %s\n', 'progress', 'elapsed', 'remaining (estimation)');
        fprintf('%08.4f%% -- %s -- %s\n', perc*100, etimev(clock,sstrt), '?');
        % pause(0.1);

    case 'done'    
        fprintf('%08.4f%% -- %s\n', 100, 'done');
    end
end
status = 0;

function S = etimev(t1,t0)
    %% etimev(t1,t0) Verbose Elapsed time.
    %   etimev(t1,t0) returns string of the days, hours, minutes, seconds that have elapsed 
    %   between vectors T1 and T0.  The two vectors must be six elements long, in
    %   the format returned by CLOCK:
    %
    %       T = [Year Month Day Hour Minute Second]
    %   OR
    %   etimev(t), t in seconds

    if (exist('t1', 'var') && exist('t0', 'var') && length(t1)>2 && length(t0)>2)
        t = etime(t1,t0);     %Time in seconds
        if t < 0
            t=-t;
        end
    elseif length(t1) == 1
        t = t1;
    else
        t = 0;
    end
    days    = floor(t/(24*60*60));
    t       = t-days*24*60*60;
    hours   = floor(t/(60*60));
    t       = t-hours*60*60;
    mins    = floor(t/60);
    t       = floor(t-mins*60);
    if days > 0
        S = [num2str(days) 'd ' num2str(hours) 'h ' num2str(mins) 'm ' num2str(t) 's'];
    elseif hours > 0
        S = [num2str(hours) 'h ' num2str(mins) 'm ' num2str(t) 's'];
    elseif mins > 0
        S = [num2str(mins) 'm ' num2str(t) 's'];
    else
        S = [num2str(t) 's'];
    end