%% DSA_run_plot
%
% DSA_run_plot is an internal script that is normally run from inside
% function DSA_plot(DSA). This script goes through the defined sets of
% variables and displays their time-domain waveforms and harmonic-domain
% spectra according to the requested HD_sources and TD_sources.
% 
% See DSA_plot(DSA) function too.
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

nb_sets = length(sets);
for iset = 1:nb_sets
    plot_vars = sets{iset};
    sz = size(plot_vars, 1);
    all_leg = {};
    
    if sz > 0
        % +++++++++++++++ HARMONIC DOMAIN +++++++++++++++
        
        if ~isempty(HD_sources)
            for ii = 1:sz % loop on the variables
                scaling_factor = 1;

                for jj = 1:size(HD_sources, 1) % loop on the HD sources
                    HD_source = DSA.HD.(HD_sources{jj, 1});
                    xx = HD_source.f;
                    yy1 = HD_source.(plot_vars{ii}) * scaling_factor;
                    abs_yy_threshold = 1e-17; %5e-3;
                    figure
                    newname = replace(plot_vars{ii},'_','');
                    sgtitle(['harmonic spectrum: '  newname ' ' HD_sources{jj, 2}], 'fontsize', 14)
                    ax1 = subplot(2, 1, 1);
                        hold on
                        grid on
                        stem(xx, abs(yy1), '.', 'linewidth', 1.2)
                        if xx(1) < xx(end)
                            xlim(xx([1, end]))
                        end
                        if max(abs(yy1)) < abs_yy_threshold
                            ylim([0, 0.5])
                        else
                            ylim([0, max(abs(yy1))*1.1])
                        end
                        ylabel('amplitude [per unit]', 'fontsize', 14)
                    ax2 = subplot(2, 1, 2);
                        hold on
                        grid on
                        stem(xx, angle(yy1)*180/pi, '.', 'linewidth', 1.2)
                        if xx(1) < xx(end)
                            xlim(xx([1, end]))
                        end
                        ylim([-190, +190])
                        xlabel('frequency [Hz]', 'fontsize', 14)
                        ylabel('phase [deg]', 'fontsize', 14)
                    linkaxes([ax1, ax2], 'x')
                end
            end
        end
        
        % +++++++++++++++ TIME DOMAIN +++++++++++++++
        
        if ~isempty(TD_sources)
            figure
            hold on
            grid on
            sgtitle('signals in time-domain', 'fontsize', 14)
            xlabel('time [s]', 'fontsize', 14)
            ylabel('quantity', 'fontsize', 14)
            ax = gca;
            
            for jj = 1:size(TD_sources, 1) % loop on the TD sources
                TD_source = DSA.TD.(TD_sources{jj, 1});
                xx = TD_source.t;
                
                ax.ColorOrderIndex = 1;
                leg = cell(1, sz);
        
                for ii = 1:sz % loop on the variables
                    scaling_factor = 1;
                    yy1 = TD_source.(plot_vars{ii}) * scaling_factor;
                    plot_density = TD_sources{jj, 4};
                    plot(xx(1:plot_density:end), yy1(1:plot_density:end), TD_sources{jj, 2}{1}, 'linewidth', 1.2);
                    newname = replace(plot_vars{ii},'_','');
                    leg{ii} = [newname ' ' TD_sources{jj, 3}];
                end
                all_leg = [all_leg leg]; %#ok<AGROW> 
            end
            
            % +++++++++++++++ LEGEND +++++++++++++++

            hleg = legend(all_leg, 'fontsize', 12, 'location', 'eastoutside');
            set(gcf, 'unit', 'inches'); set(hleg, 'unit', 'inches');
            figure_size = get(gcf, 'position'); legend_size = get(hleg, 'position');
            figure_size(3) = figure_size(3) + legend_size(3);
            set(gcf, 'position', figure_size)
            
            % +++++++++++++++ Y-AXIS automatic expansion +++++++++++++++
        
            yl = ylim;
            if abs(yl(1)-yl(2)) < 5e-4
               ylim(yl(1) + [-0.5 +0.5])
            end
        end
    else
        warning('DSA_plot: skips empty set')
    end
end