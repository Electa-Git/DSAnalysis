function DSA = DSA_run_FRESP(DSA, TFid)
%% DSA = DSA_run_FRESP(DSA, TFid)
% 
% This function calculate the frequency response of the complete lifted
% system DSA.LTI.HSS_SYS, by calling freqresp on the HSS system over a 
% vector of frequency values f_values.
%
% The frequency response is calculated between input 'ivarname' and output
% 'overname'.
% The SISO frequency response between this input-output pair is extracted
% from the HTF in the sense of harmonic linearisation (HLIN) between input
% harmonic channel h_in and output harmonic channel h_out = h_lin + h_shift.
%
% sign factor is 1 or -1, to switch between load and generator convention
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

    m   = DSA.info.m;    % unlifted input  size
    p   = DSA.info.p;    % unlifted output size
    h_t = DSA.FLIFT.h_t; % lifted truncation rank
    
    % variable and ONE harmonic shift.
    TF_info     = DSA.FRESP.(TFid);
    f_values    = TF_info.f_values;
    ivarname    = TF_info.ivarname;
    ovarname    = TF_info.ovarname;
    h_in        = TF_info.h_in;
    h_shift     = TF_info.h_shift;
    sign_factor = TF_info.sign_factor;
    sys         = DSA.LTI.HSS_SYS;
    
    % -- extract frequency response of FULL system
    f_resp_HTF = freqresp(sys, f_values, 'Hz');

    % -- get index of input and output variable
    in_idx_unlifted  = DSA.vars.lidx.smsi_inputs.(ivarname);
    out_idx_unlifted = DSA.vars.lidx.smsi_outputs.(ovarname);
    
    % -- mapping from variable index and harmonic channel to frequency-lifted index:
    h_out = h_in+h_shift; % output harmonic channel
    in_idx  = (h_in  + h_t).*m + in_idx_unlifted ;
    out_idx = (h_out + h_t).*p + out_idx_unlifted;
    
    % -- extract HLIN frequency response of the requested i-o pair
    f_resp_HLIN = sign_factor*squeeze(f_resp_HTF(out_idx, in_idx, :));
    
    % -- output results
    DSA.FRESP.(TFid).HTF.f_resp  = f_resp_HTF;
    DSA.FRESP.(TFid).HLIN.f_resp = f_resp_HLIN;

    %% -- plot magnitude and phase over frequency range
    % useful functions: rad2deg, unwrap, wrapTo360, wrapTo180
    
    ivarname_str = replace(ivarname, '_', '');
    ovarname_str = replace(ovarname, '_', '');

    figure
    sgtitle(sprintf('HLIN freqresp: input %s to output %s', ivarname_str, ovarname_str))
    ax1 = subplot(2, 1, 1);
        hold on
        grid on
        xx = f_values;
        yy = 20*log10(abs(f_resp_HLIN));
        plot(xx, yy, 'linewidth', 1.5);
        ylims = ylim;
        set(gca, 'XScale', 'log')
        ylabel('magnitude [dB]', 'Fontsize', 14)
        xlim([min(xx([1, end])), max(xx([1, end]))])
        ylim([ylims(1)-5, ylims(2)+5])

    ax2 = subplot(2, 1, 2);
        hold on
        grid on
        xx = f_values;
        yy = rad2deg(angle(f_resp_HLIN));
        plot(xx, yy, 'linewidth', 1.5)
        ylims = ylim;
        set(gca, 'XScale', 'log')
        xlabel('frequency [Hz]', 'Fontsize', 14)
        ylabel('phase [deg]', 'Fontsize', 14)
        xlim([min(xx([1, end])), max(xx([1, end]))])
        ylim([ylims(1)-5, ylims(2)+5])
        yline([-90  90], 'color','#D95319')
    linkaxes([ax1, ax2],'x')

end
