function DSA = DSA_run_FRESP_3PH(DSA, TFid)
%% DSA = DSA_run_FRESP_3PH(DSA, TFid)
% 
% This function calculate the frequency response of the complete lifted
% system DSA.LTI.HSS_SYS, by calling freqresp on the HSS system over a 
% vector of frequency values f_values.
%
% This function is a variation of DSA_run_FRESP specifically dedicated to
% three-phase transfer functions.
% The frequency response is calculated between inputs 'ivarname' and
% outputs 'overname'.
% The MIMO frequency response between this input-output pair is extracted
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
    in_idx_a = DSA.vars.lidx.smsi_inputs.([ivarname '_a']);
    in_idx_b = DSA.vars.lidx.smsi_inputs.([ivarname '_b']);
    in_idx_c = DSA.vars.lidx.smsi_inputs.([ivarname '_c']);
    in_idx_unlifted = [in_idx_a, in_idx_b, in_idx_c];

    out_idx_a = DSA.vars.lidx.smsi_outputs.([ovarname '_a']);
    out_idx_b = DSA.vars.lidx.smsi_outputs.([ovarname '_b']);
    out_idx_c = DSA.vars.lidx.smsi_outputs.([ovarname '_c']);
    out_idx_unlifted = [out_idx_a, out_idx_b, out_idx_c];

    % -- extract frequency response of the requested i-o pair
    
    % -- mapping from variable index and harmonic channel to frequency-lifted index:
    h_out = h_in+h_shift; % output harmonic channel
    in_idx  = (h_in  + h_t).*m + in_idx_unlifted ;
    out_idx = (h_out + h_t).*p + out_idx_unlifted;
    
    % -- extract HLIN frequency response of the requested i-o pair
    f_resp_HLIN_abc = sign_factor*squeeze(f_resp_HTF(out_idx, in_idx, :)); % sign factor is 1 or -1, so that TF definition agrees with gen or load sign conventions.
    
    f_resp_HLIN_a2a = squeeze(f_resp_HLIN_abc(1,1,:));
    f_resp_HLIN_a2b = squeeze(f_resp_HLIN_abc(2,1,:));
    f_resp_HLIN_a2c = squeeze(f_resp_HLIN_abc(3,1,:));
    f_resp_HLIN_b2a = squeeze(f_resp_HLIN_abc(1,2,:));
    f_resp_HLIN_b2b = squeeze(f_resp_HLIN_abc(2,2,:));
    f_resp_HLIN_b2c = squeeze(f_resp_HLIN_abc(3,2,:));
    f_resp_HLIN_c2a = squeeze(f_resp_HLIN_abc(1,3,:));
    f_resp_HLIN_c2b = squeeze(f_resp_HLIN_abc(2,3,:));
    f_resp_HLIN_c2c = squeeze(f_resp_HLIN_abc(3,3,:));
    
    % -- output results
    DSA.FRESP.(TFid).HTF.fresp         = f_resp_HTF;
    DSA.FRESP.(TFid).HLIN.MIMO.abc2abc = f_resp_HLIN_abc;
    DSA.FRESP.(TFid).HLIN.SISO.a2a     = f_resp_HLIN_a2a;
    DSA.FRESP.(TFid).HLIN.SISO.a2b     = f_resp_HLIN_a2b;
    DSA.FRESP.(TFid).HLIN.SISO.a2c     = f_resp_HLIN_a2c;
    DSA.FRESP.(TFid).HLIN.SISO.b2a     = f_resp_HLIN_b2a;
    DSA.FRESP.(TFid).HLIN.SISO.b2b     = f_resp_HLIN_b2b;
    DSA.FRESP.(TFid).HLIN.SISO.b2c     = f_resp_HLIN_b2c;
    DSA.FRESP.(TFid).HLIN.SISO.c2a     = f_resp_HLIN_c2a;
    DSA.FRESP.(TFid).HLIN.SISO.c2b     = f_resp_HLIN_c2b;
    DSA.FRESP.(TFid).HLIN.SISO.c2c     = f_resp_HLIN_c2c;

    %% -- transformation to sequences
    
    f_resp_HLIN_zpn = zeros(size(f_resp_HLIN_abc));
    a = exp(1j*2*pi/3);
    T = [1, 1, 1;
         1, a^2, a;
         1, a, a^2]; % Fortescue
    
    % applying Fortescue at each frequency:
    for idx = 1:size(f_resp_HLIN_abc, 3)
        f_resp_HLIN_zpn(:,:,idx) = T\f_resp_HLIN_abc(:,:,idx)*T;
    end
    
    % inputs are columns, second index
    % outputs are rows, first index
    f_resp_HLIN_z2z = squeeze(f_resp_HLIN_zpn(1,1,:));
    f_resp_HLIN_z2p = squeeze(f_resp_HLIN_zpn(2,1,:));
    f_resp_HLIN_z2n = squeeze(f_resp_HLIN_zpn(3,1,:));
    f_resp_HLIN_p2z = squeeze(f_resp_HLIN_zpn(1,2,:));
    f_resp_HLIN_p2p = squeeze(f_resp_HLIN_zpn(2,2,:));
    f_resp_HLIN_p2n = squeeze(f_resp_HLIN_zpn(3,2,:));
    f_resp_HLIN_n2z = squeeze(f_resp_HLIN_zpn(1,3,:));
    f_resp_HLIN_n2p = squeeze(f_resp_HLIN_zpn(2,3,:));
    f_resp_HLIN_n2n = squeeze(f_resp_HLIN_zpn(3,3,:));
    
    % -- output results
    DSA.FRESP.(TFid).HLIN.MIMO.zpn2zpn = f_resp_HLIN_zpn;
    DSA.FRESP.(TFid).HLIN.SISO.z2z     = f_resp_HLIN_z2z;
    DSA.FRESP.(TFid).HLIN.SISO.z2p     = f_resp_HLIN_z2p;
    DSA.FRESP.(TFid).HLIN.SISO.z2n     = f_resp_HLIN_z2n;
    DSA.FRESP.(TFid).HLIN.SISO.p2z     = f_resp_HLIN_p2z;
    DSA.FRESP.(TFid).HLIN.SISO.p2p     = f_resp_HLIN_p2p;
    DSA.FRESP.(TFid).HLIN.SISO.p2n     = f_resp_HLIN_p2n;
    DSA.FRESP.(TFid).HLIN.SISO.n2z     = f_resp_HLIN_n2z;
    DSA.FRESP.(TFid).HLIN.SISO.n2p     = f_resp_HLIN_n2p;
    DSA.FRESP.(TFid).HLIN.SISO.n2n     = f_resp_HLIN_n2n;

    %% -- plot magnitude and phase over frequency range: all not involving zero sequence
    % useful functions: rad2deg, unwrap, wrapTo360, wrapTo180

    ivarname_str = replace(ivarname, '_', '');
    ovarname_str = replace(ovarname, '_', '');
    hfig = figure;
    sgtitle(sprintf('TF: input %s to output %s', ivarname_str, ovarname_str))
    
    leg = {};
    ax1 = subplot(2, 1, 1);
        hold on
        grid on
        xx = f_values;
        if strcmp(get(groot, 'defaultTextInterpreter'), 'latex')
            locstr = '$';
        else
            locstr = '';
        end
        plot(xx, 20*log10(abs(f_resp_HLIN_p2p)), '-' , 'linewidth', 1.2); leg = [leg, [locstr 'p \rightarrow p' locstr]];
        plot(xx, 20*log10(abs(f_resp_HLIN_n2n)), '-' , 'linewidth', 1.2); leg = [leg, [locstr 'n \rightarrow n' locstr]];
        plot(xx, 20*log10(abs(f_resp_HLIN_z2z)), '--', 'linewidth', 0.5); leg = [leg, [locstr 'z \rightarrow z' locstr]];
        plot(xx, 20*log10(abs(f_resp_HLIN_z2p)), '--', 'linewidth', 0.5); leg = [leg, [locstr 'z \rightarrow p' locstr]];
        plot(xx, 20*log10(abs(f_resp_HLIN_z2n)), '--', 'linewidth', 0.5); leg = [leg, [locstr 'z \rightarrow n' locstr]];
        plot(xx, 20*log10(abs(f_resp_HLIN_p2z)), '--', 'linewidth', 0.5); leg = [leg, [locstr 'p \rightarrow z' locstr]];
        plot(xx, 20*log10(abs(f_resp_HLIN_n2z)), '--', 'linewidth', 0.5); leg = [leg, [locstr 'n \rightarrow z' locstr]];
        plot(xx, 20*log10(abs(f_resp_HLIN_p2n)), '--', 'linewidth', 0.5); leg = [leg, [locstr 'p \rightarrow n' locstr]];
        plot(xx, 20*log10(abs(f_resp_HLIN_n2p)), '--', 'linewidth', 0.5); leg = [leg, [locstr 'n \rightarrow p' locstr]];
        set(gca, 'XScale', 'log')
        ylabel('magnitude [dB]', 'Fontsize', 14)
        legend(leg, 'Fontsize', 14, 'Location', 'Eastoutside')
        ax1.Position = [0.1300    0.5635    0.5750    0.3293];

    hfig.Position = [680   558   1.5*560   420];
    ax2 = subplot(2, 1, 2);
        hold on
        grid on
        xx = f_values;
        plot(xx, rad2deg(angle(f_resp_HLIN_p2p)), '-' , 'linewidth', 1.2)
        plot(xx, rad2deg(angle(f_resp_HLIN_n2n)), '-' , 'linewidth', 1.2)
        plot(xx, rad2deg(angle(f_resp_HLIN_z2z)), '--', 'linewidth', 0.5)
        plot(xx, rad2deg(angle(f_resp_HLIN_z2p)), '--', 'linewidth', 0.5)
        plot(xx, rad2deg(angle(f_resp_HLIN_z2n)), '--', 'linewidth', 0.5)
        plot(xx, rad2deg(angle(f_resp_HLIN_p2z)), '--', 'linewidth', 0.5)
        plot(xx, rad2deg(angle(f_resp_HLIN_n2z)), '--', 'linewidth', 0.5)
        plot(xx, rad2deg(angle(f_resp_HLIN_p2n)), '--', 'linewidth', 0.5)
        plot(xx, rad2deg(angle(f_resp_HLIN_n2p)), '--', 'linewidth', 0.5)
        yline([-90  90], 'color','#D95319')
        set(gca, 'XScale', 'log')
        xlabel('frequency [Hz]', 'Fontsize', 14)
        ylabel('phase [deg]', 'Fontsize', 14)
        ax2.Position = [0.1300    0.1133    0.5750    0.3311];
        
    linkaxes([ax1, ax2],'x')
end
