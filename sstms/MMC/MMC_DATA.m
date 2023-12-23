function DSA = MMC_DATA(DSA)
    % ---------- GENERAL NETWORK DATA ---------- %
    % source: Julian Freytes PhD Thesis
    
    f_1                 = 50;       % [Hz], fundamental frequency
    omega_1             = 2*pi*f_1; % [rad/s]
    S_nom_MVA           = 1000;     % [MVA]
    U_a_nom_LL_kVrms    = 320;      % [kVrms]
    V_d_nom_kV          = 640;      % [kV]
    
    % ---------- BASE VALUES CALCULATION ---------- %
    
    S_base   = S_nom_MVA*1e6;                   % VA
    E_base   = 2/3 * S_base;                    % VA
    V_a_base = U_a_nom_LL_kVrms*1e3*sqrt(2/3);  % phase peak V
    I_a_base = 2/3 * S_base / V_a_base;         % phase peak A
    Z_a_base = V_a_base / I_a_base;             % Ohm
    V_d_nom  = V_d_nom_kV * 1e3;                % V
    I_d_nom  = S_base / V_d_nom;                % A
    
    % These base values are used for both AC and DC sides
    base.S_base   = S_base;
    base.E_base   = E_base;
    base.V_a_base = V_a_base;
    base.I_a_base = I_a_base;
    base.Z_a_base = Z_a_base;
    
    % ---------- TRANSFORMATION RATIOS ---------- %
    
    % none
    
    % ---------- CIRCUIT PARAMETERS ---------- %
    % source: Julian Freytes PhD Thesis
    
    X_g_pu = 0.18;
    R_g_pu = 0.005;
    R_a_pu = 0.010;
    C_s    = 13e-3; % [F] 
    N_s    = 400;
    C_a    = C_s/N_s;
    L_a    = 48e-3; % [F]
    
    L_g    = X_g_pu * Z_a_base / omega_1;
    R_g    = R_g_pu * Z_a_base;
    R_a    = R_a_pu * Z_a_base;
    
    L_e    = L_a/2 + L_g;
    R_e    = R_a/2 + R_g;
    
    param.R_a = R_a / Z_a_base;
    param.L_a = L_a / Z_a_base;
    param.C_a = C_a * Z_a_base;
    param.R_g = R_g / Z_a_base;
    param.L_g = L_g / Z_a_base;
    param.R_e = R_e / Z_a_base;
    param.L_e = L_e / Z_a_base;
    
    param.omega_1 = omega_1 / 1;
    
    % ---------- CONTROL PARAMETERS ---------- %
    
    % AC: MMC AC current controller
    % CC: MMC circulating current controller
    % DC: MMC DC voltage controller
    % AE: MMC arm-energy balancing controller
    % NF: MMC notch filter
    % PLL: phase locked loop
    
    BW_AC  = 2*pi*150;
    BW_CC  = 2*pi*150;
    BW_PQ  = 2*pi*5;
    BW_PLL = 2*pi*10;
    BW_AE  = 2*pi*10;
    BW_NF  = 2*pi*10;
    
    AC_Kp  = BW_AC * L_e;
    AC_Kr  = BW_AC * R_e;
    CC_Kp  = BW_CC * L_a;
    CC_Kr  = BW_CC * R_a;
    AE_Kp  = 2*sqrt(2)/3*BW_AE;
    AE_Ki  = 1/9*BW_AE;
    PLL_Kp = 2*sqrt(2)/3*BW_PLL;
    PLL_Ki = 1/9*BW_PLL^2;
    PQ_Ki  = 2/3/V_a_base * BW_PQ;
    PQ_Kp  = 2/3/V_a_base * 0.05; % 0.05 << 1
    
    param.AC_Kp    = AC_Kp / Z_a_base;
    param.AC_Kr    = AC_Kr / Z_a_base;
    param.CC_Kp    = CC_Kp / Z_a_base;
    param.CC_Kr    = CC_Kr / Z_a_base;
    param.PQ_Kp    = PQ_Kp * V_a_base;
    param.PQ_Ki    = PQ_Ki * V_a_base;
    param.AE_Kp    = AE_Kp / 1;
    param.AE_Ki    = AE_Ki / 1;
    param.PLL_Kp   = PLL_Kp / 1;
    param.PLL_Ki   = PLL_Ki / 1;
    param.omega_NF = BW_NF;
    
    % ---------- EXTERNAL HD-SOURCES ---------- %
    
    % (none)

    % ---------- INPUT PARAMETERS ---------- %
    
    P_g_ref      =  0.8 * S_base;
    Q_g_ref      = -0.3 * S_base;
    P_g_ref_step = 0.85 * S_base;
    Q_g_ref_step = -0.3 * S_base;
    W_Sbar_ref   = C_a * V_d_nom^2;
    W_Dbar_ref   = 0;
    V_g_ph_pk    = V_a_base;
    I_g_ph_pk    = I_a_base;
    I_a_ph_pk    = I_a_base;
    V_d_pk       = V_d_nom;
    I_d_pk       = I_d_nom;
    
    param.V_g_ph_pk    = V_g_ph_pk  / V_a_base;
    param.I_g_ph_pk    = I_g_ph_pk  / I_a_base;
    param.I_a_ph_pk    = I_a_ph_pk  / I_a_base;
    param.V_d_pk       = V_d_pk     / V_a_base;
    param.I_d_pk       = I_d_pk     / I_a_base;
    param.S_nom        = S_base     / E_base;
    param.P_g_ref      = P_g_ref    / E_base;
    param.Q_g_ref      = Q_g_ref    / E_base;
    param.P_g_ref_step = P_g_ref_step / E_base;
    param.Q_g_ref_step = Q_g_ref_step / E_base;
    param.W_Sbar_ref   = W_Sbar_ref / E_base;
    param.W_Dbar_ref   = W_Dbar_ref / E_base;
    
    % ---------- DELAY PARAMETERS ---------- %
    % setting a non-zero delay value only impacts the systems implemented
    % in such a way that they account for delays. The non-delayed systems
    % disregard the DSA.delay field.
    
    T_d_us = 150;
    T_d = T_d_us*1e-6;
    
    delay.n_u_a_ref = T_d;
    delay.n_u_b_ref = T_d;
    delay.n_u_c_ref = T_d;
    delay.n_l_a_ref = T_d;
    delay.n_l_b_ref = T_d;
    delay.n_l_c_ref = T_d;
    
    delay.n_u_ref = T_d;
    delay.n_l_ref = T_d;

    % ---------- PARAMETERS STRUCTURES ---------- %
    
    DSA.data.f_1   = f_1;
    DSA.data.delay = delay;
    DSA.data.base  = base;
    DSA.data.param = param;
end