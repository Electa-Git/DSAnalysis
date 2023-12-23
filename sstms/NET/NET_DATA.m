function DSA = NET_DATA(DSA)
    % ---------- GENERAL NETWORK DATA ---------- %
    % source for the MMC and HVDC: Julian Freytes PhD Thesis (except for G_a)
    % source for the medium-voltage TLCs: Yiran Jing PhD thesis
    %
    % "primary subscripts":
    %
    % a: MMCs AC terminals, on the converter side of the transformer
    % d: refers to direct currents or voltages, not the direct component of
    %    dq systems!
    % g: MMCs AC terminals, on the grid side of the transformer
    % t: related to the TLCs transformers
    % v: related to MMCs or to the HVDC link in general
    % w: related to TLCs or to the collection grid in general
    % x: MV cables of the collection grid
    %
    % "secondary subscripts":
    % Some primary subscripts are used as secondary subscripts too.
    % 
    % a: phase a
    % b: phase b
    % c: phase c
    % al: alpha component
    % be: beta component
    % d: direct component of a dq system
    % q: quadrature component of a dq system
    % flt: filtered
    % quad: in quadrature (SOGI-PNSE)
    % ps: positive sequence
    % ns: negative sequence
    % zs: zero sequence
    
    f_1                = 50;       % [Hz]    fundamental frequency
    omega_1            = 2*pi*f_1; % [rad/s] 
    S_v_nom_MVA        = 1000;     % [MVA]   whole system rated apparent power
    S_w_1_nom_MVA      = 500;      % [MVA]   TLC1 rated apparent power
    S_w_2_nom_MVA      = 400;      % [MVA]   TLC2 rated apparent power
    S_w_3_nom_MVA      = 100;      % [MVA]   TLC3 rated apparent power
    U_g_1_nom_LL_kVrms = 380;      % [kVrms] onshore grid voltage
    U_g_2_nom_LL_kVrms = 200;      % [kVrms] collection grid voltage
    U_a_nom_LL_kVrms   = 320;      % [kVrms] MMC AC terminal voltage
    U_w_nom_LL_kVrms   = 33;       % [kVrms] TLC AC terminal voltage
    V_d_nom_kV         = 640;      % [kV]    MMC DC terminal voltage
    V_d_w_nom_kV       = 65;       % [kV]    TLC DC terminal voltage
    
    % ---------- BASE VALUES CALCULATION ---------- %
    
    S_base     = S_v_nom_MVA*1e6;                 % VA
    E_base     = 2/3 * S_base;                    % VA
    V_a_base   = U_a_nom_LL_kVrms*1e3*sqrt(2/3);  % phase peak V
    I_a_base   = E_base / V_a_base;               % phase peak A
    Z_a_base   = V_a_base / I_a_base;             % Ohm
    V_g_1_base = U_g_1_nom_LL_kVrms*1e3*sqrt(2/3); 
    I_g_1_base = E_base / V_g_1_base;          
    Z_g_1_base = V_g_1_base / I_g_1_base;              
    V_g_2_base = U_g_2_nom_LL_kVrms*1e3*sqrt(2/3);
    I_g_2_base = E_base / V_g_2_base;
    Z_g_2_base = V_g_2_base / I_g_2_base;
    V_w_base   = U_w_nom_LL_kVrms*1e3*sqrt(2/3);
    I_w_base   = E_base / V_w_base;
    Z_w_base   = V_w_base / I_w_base;
    
    S_w_1_base = S_w_1_nom_MVA*1e6;
    S_w_2_base = S_w_2_nom_MVA*1e6;
    S_w_3_base = S_w_3_nom_MVA*1e6;
    I_w_1_base = 2/3 * S_w_1_base / V_w_base;
    I_w_2_base = 2/3 * S_w_2_base / V_w_base;
    I_w_3_base = 2/3 * S_w_3_base / V_w_base;
    Z_w_1_base = V_w_base / I_w_1_base;
    Z_w_2_base = V_w_base / I_w_2_base;
    Z_w_3_base = V_w_base / I_w_3_base;
    
    V_d_nom   = V_d_nom_kV * 1e3;
    I_d_nom   = S_base / V_d_nom;
    V_d_w_nom = V_d_w_nom_kV * 1e3;
    
    base.S_base     = S_base;
    base.E_base     = E_base;
    base.V_a_base   = V_a_base;
    base.I_a_base   = I_a_base;
    base.Z_a_base   = Z_a_base;
    base.V_g_1_base = V_g_1_base;
    base.I_g_1_base = I_g_1_base;
    base.Z_g_1_base = Z_g_1_base;
    base.V_g_2_base = V_g_2_base;
    base.I_g_2_base = I_g_2_base;
    base.Z_g_2_base = Z_g_2_base;
    base.V_w_base   = V_w_base;
    base.I_w_base   = I_w_base;
    base.Z_w_base   = Z_w_base;
    base.V_d_nom    = V_d_nom;
    base.I_d_nom    = I_d_nom;
    
    param.omega_1 = omega_1;
    
    % ---------- TRANSFORMATION RATIOS ---------- %
    
    % Converter terminal voltage always taken as denominator, for
    % consistency. Consequently, the transformation ratio may be greater or
    % lower than 1.
    k_MMC1  = V_g_1_base/V_a_base;
    k_MMC2  = V_g_2_base/V_a_base;
    k_TLC   = V_g_2_base/V_w_base;
    
    param.k_MMC1 = k_MMC1/k_MMC1; % unity in the scaled system
    param.k_MMC2 = k_MMC2/k_MMC2; % unity in the scaled system
    param.k_TLC  = k_TLC/k_TLC;   % unity in the scaled system
    
    % ---------- CIRCUIT PARAMETERS ---------- %
    
    % --- MMCs and their transformers ---
    
    X_g_pu = 0.18;
    R_g_pu = 0.005;
    R_a_pu = 0.01;

    N_s_1 = 400;
    C_s_1 = 13e-3; % [F]
    C_a_1 = C_s_1/N_s_1;
    R_a_1 = R_a_pu * Z_a_base;
    G_a_1 = 1e-9;  % [S]
    L_a_1 = 48e-3; % [F]
    
    N_s_2 = 400;
    C_s_2 = 13e-3; % [F]
    C_a_2 = C_s_2/N_s_2;
    R_a_2 = R_a_pu * Z_a_base;
    G_a_2 = 1e-9;  % [S]
    L_a_2 = 48e-3; % [F]
    
    L_g_1 = X_g_pu * Z_a_base / omega_1;
    R_g_1 = R_g_pu * Z_a_base;
    L_e_1 = L_a_1/2 + L_g_1;
    R_e_1 = R_a_1/2 + R_g_1;
    
    L_g_2 = X_g_pu * Z_a_base / omega_1;
    R_g_2 = R_g_pu * Z_a_base;
    L_e_2 = L_a_2/2 + L_g_2;
    R_e_2 = R_a_2/2 + R_g_2;
    
    % --- TLC transformer and filter ---
    Xl_f_pu = 0.2;
    Xc_f_pu = 1/0.15;
    R_f_pu  = 0.01;
    
    L_f_1 = Xl_f_pu * Z_w_1_base / omega_1;
    L_f_2 = Xl_f_pu * Z_w_2_base / omega_1;
    L_f_3 = Xl_f_pu * Z_w_3_base / omega_1;
    R_f_1 = R_f_pu * Z_w_1_base;
    R_f_2 = R_f_pu * Z_w_2_base;
    R_f_3 = R_f_pu * Z_w_3_base;
    C_f_1 = 1/(Xc_f_pu * Z_w_1_base) / omega_1;
    C_f_2 = 1/(Xc_f_pu * Z_w_2_base) / omega_1;
    C_f_3 = 1/(Xc_f_pu * Z_w_3_base) / omega_1;
    
    X_t_pu = 0.1;
    R_t_pu = 0.01;
    
    L_t_1 = X_t_pu * Z_w_1_base / omega_1;
    L_t_2 = X_t_pu * Z_w_2_base / omega_1;
    L_t_3 = X_t_pu * Z_w_3_base / omega_1;
    R_t_1 = R_t_pu * Z_w_1_base;
    R_t_2 = R_t_pu * Z_w_2_base;
    R_t_3 = R_t_pu * Z_w_3_base;
    
    % --- cables ---
    cable_0_km = 100; % [km]
    cable_1_km = 10;  % [km]
    cable_2_km = 5;   % [km]
    cable_3_km = 5;   % [km]
    
    R_d_pkm   = 15e-3;    % [ohm/km]
    L_d_pkm   = 0.3e-3;   % [H/km]
    C_d_pkm   = 0.12e-6;  % [F/km]
    G_d_pkm   = 0;        % [S/km]
    R_x_1_pkm = 15e-3;    % [ohm/km]
    R_x_2_pkm = 16.5e-3;  % [ohm/km]
    R_x_3_pkm = 16.5e-3;  % [ohm/km]
    L_x_1_pkm = 0.3e-3;   % [H/km]
    L_x_2_pkm = 0.33e-3;  % [H/km]
    L_x_3_pkm = 0.33e-3;  % [H/km]
    C_x_1_pkm = 3*0.12e-6;  % [F/km]
    C_x_2_pkm = 3*0.11e-6;  % [F/km]
    C_x_3_pkm = 3*0.11e-6;  % [F/km]
    G_x = 1e-8; % [S], small fixed value to avoid poles on the imaginary axis
    
    R_d   = cable_0_km * R_d_pkm;
    L_d   = cable_0_km * L_d_pkm;
    C_d   = cable_0_km * C_d_pkm; % pole to neutral capacitance
    G_d   = cable_0_km * G_d_pkm;
    R_x_1 = cable_1_km * R_x_1_pkm;
    L_x_1 = cable_1_km * L_x_1_pkm;
    C_x_1 = cable_1_km * C_x_1_pkm; % phase to neutral capacitance
    R_x_2 = cable_2_km * R_x_2_pkm;
    L_x_2 = cable_2_km * L_x_2_pkm;
    C_x_2 = cable_2_km * C_x_2_pkm; % phase to neutral capacitance
    R_x_3 = cable_3_km * R_x_3_pkm;
    L_x_3 = cable_3_km * L_x_3_pkm;
    C_x_3 = cable_3_km * C_x_3_pkm; % phase to neutral capacitance
    % f_res_1 = 1/sqrt(L_x_1*C_x_1)
    % f_res_2 = 1/sqrt(L_x_2*C_x_2)
    % f_res_3 = 1/sqrt(L_x_3*C_x_3)
    
    param.R_a_1 = R_a_1 / Z_a_base;
    param.L_a_1 = L_a_1 / Z_a_base;
    param.C_a_1 = C_a_1 * Z_a_base;
    param.G_a_1 = G_a_1 * Z_a_base;
    param.R_a_2 = R_a_2 / Z_a_base;
    param.L_a_2 = L_a_2 / Z_a_base;
    param.C_a_2 = C_a_2 * Z_a_base;
    param.G_a_2 = G_a_2 * Z_a_base;
    param.R_g_1 = R_g_1 / Z_a_base;
    param.L_g_1 = L_g_1 / Z_a_base;
    param.R_e_1 = R_e_1 / Z_a_base;
    param.L_e_1 = L_e_1 / Z_a_base;
    param.R_g_2 = R_g_2 / Z_a_base;
    param.L_g_2 = L_g_2 / Z_a_base;
    param.R_e_2 = R_e_2 / Z_a_base;
    param.L_e_2 = L_e_2 / Z_a_base;
    param.R_d   = R_d   / Z_a_base;
    param.L_d   = L_d   / Z_a_base;
    param.C_d   = C_d   * Z_a_base;
    param.G_d   = G_d   * Z_a_base;

    param.R_f_1 = R_f_1 / Z_w_base;
    param.R_f_2 = R_f_2 / Z_w_base;
    param.R_f_3 = R_f_3 / Z_w_base;
    param.L_f_1 = L_f_1 / Z_w_base;
    param.L_f_2 = L_f_2 / Z_w_base;
    param.L_f_3 = L_f_3 / Z_w_base;
    param.C_f_1 = C_f_1 * Z_w_base;
    param.C_f_2 = C_f_2 * Z_w_base;
    param.C_f_3 = C_f_3 * Z_w_base;
    param.R_t_1 = R_t_1 / Z_w_base;
    param.R_t_2 = R_t_2 / Z_w_base;
    param.R_t_3 = R_t_3 / Z_w_base;
    param.L_t_1 = L_t_1 / Z_w_base;
    param.L_t_2 = L_t_2 / Z_w_base;
    param.L_t_3 = L_t_3 / Z_w_base;
    
    param.R_x_1 = R_x_1 / Z_g_2_base;
    param.R_x_2 = R_x_2 / Z_g_2_base;
    param.R_x_3 = R_x_3 / Z_g_2_base;
    param.L_x_1 = L_x_1 / Z_g_2_base;
    param.L_x_2 = L_x_2 / Z_g_2_base;
    param.L_x_3 = L_x_3 / Z_g_2_base;
    param.C_x_1 = C_x_1 * Z_g_2_base;
    param.C_x_2 = C_x_2 * Z_g_2_base;
    param.C_x_3 = C_x_3 * Z_g_2_base;
    param.G_x   = G_x   * Z_g_2_base;
    
    % ---------- EXTERNAL SYSTEM ---------- %
    
    % HVDC cable (fitted state-space):
    
    param.rr1 = 2.67973062e+02 / Z_a_base;
    param.rr2 = 2.16953994e+01 / Z_a_base;
    param.rr3 = 1.90508581e+00 / Z_a_base;
    param.rr4 = 1.58545014e+00 / Z_a_base;
    param.rr5 = 1.76733416e-01 / Z_a_base;
    param.rr6 =-2.14558031e-01 / Z_a_base;
    param.ll1 = 1.33947266e-02 / Z_a_base;
    param.ll2 = 7.19239114e-03 / Z_a_base;
    param.ll3 = 5.18656381e-03 / Z_a_base;
    param.ll4 = 4.24290724e-02 / Z_a_base;
    param.ll5 = 5.42725024e-02 / Z_a_base;
    param.ll6 =-8.95099527e+03 / Z_a_base;
    param.cc  = 3.99449696e-06 * Z_a_base;
    param.gg  = (6.00000002e+05)^-1 * Z_a_base;

    % ---------- CONTROL PARAMETERS ---------- %
    
    % BW:  bandwidth [rad/s]
    % AC:  alternating current controller
    % MAC: alternating current controller of the MMC
    % TAC: alternating current controller of the TLC
    % AV:  alternating voltage controller
    % CC:  circulating current controller
    % DC:  direct voltage controller
    % PQ:  active and reactive power controller
    % MPQ: active and reactive power controller of the MMC
    % TPQ: active and reactive power controller of the TLC
    % PLL: phase locked loop
    
    BW_PLL   = 2*pi*10;
    BW_MPQ   = 2*pi*5;
    BW_MAC_1 = 2*pi*150;
    BW_MAC_2 = 2*pi*300;
    BW_MCC_1 = 2*pi*150;
    BW_MCC_2 = 2*pi*150;
    BW_MDV   = 2*pi*5;
    % BW_MDV   = 2*pi*65; warning('using unstable MDV bandwidth')
    BW_MAV   = 2*pi*200; % original
    BW_TPQ   = 2*pi*5;
    BW_TAC   = 2*pi*150;
    
    %warning('considering unstable MMC case')
    MAC_Kp_1 = BW_MAC_1 * L_e_1; % collocation paper: *-1 % / 15;
    MAC_Kr_1 = BW_MAC_1 * R_e_1; % collocation paper: *-1 % / 15;
    MAC_Kp_2 = BW_MAC_2 * L_e_2;
    MAC_Kr_2 = BW_MAC_2 * R_e_2;
    
    TAC_Kp_1 = BW_TAC * L_f_1;
    TAC_Kr_1 = BW_TAC * R_f_1;
    TAC_Kp_2 = BW_TAC * L_f_2;
    TAC_Kr_2 = BW_TAC * R_f_2;
    TAC_Kp_3 = BW_TAC * L_f_3;
    TAC_Kr_3 = BW_TAC * R_f_3;
    
    MCC_Kp_1 = BW_MCC_1 * L_a_1;
    MCC_Kr_1 = BW_MCC_1 * R_a_1;
    MCC_Kp_2 = BW_MCC_2 * L_a_2;
    MCC_Kr_2 = BW_MCC_2 * R_a_2;
    
    PLL_Kp       = 2*sqrt(2)/3*BW_PLL;
    PLL_Ki       = 1/9*BW_PLL^2;
    
    MDV_Kp_1     = 2*sqrt(2)/3*BW_MDV;
    MDV_Ki_1     = 1/9*BW_MDV^2;
    MDV_C_d_1_eq = 6*C_a_1 + (C_d/2)/2; % Equation (3.235) in Kamran's book
    
    MAV_C_g_2_eq = 1/2 * (C_x_1 + C_x_2 + C_x_3);
    MAV_Kp       = BW_MAV * MAV_C_g_2_eq;
    MAV_Kr       = 0.01; % should be ≈ 0
    
    MPQ_Ki_1 = 2/3/V_a_base * BW_MPQ;
    MPQ_Kp_1 = 2/3/V_a_base * 0.05; % 0.05 << 1
    MPQ_Ki_2 = 2/3/V_a_base * BW_MPQ;
    MPQ_Kp_2 = 2/3/V_a_base * 0.05; % 0.05 << 1
    TPQ_Ki   = 2/3/V_w_base * BW_TPQ;
    TPQ_Kp   = 2/3/V_w_base * 0.05; % 0.05 << 1
    
    param.MAV_Kp   = MAV_Kp * Z_g_2_base;
    param.MAV_Kr   = MAV_Kr * Z_g_2_base;
    param.MAC_Kp_1 = MAC_Kp_1 / Z_a_base;
    param.MAC_Kr_1 = MAC_Kr_1 / Z_a_base;
    param.MAC_Kp_2 = MAC_Kp_2 / Z_a_base;
    param.MAC_Kr_2 = MAC_Kr_2 / Z_a_base;
    param.MCC_Kp_1 = MCC_Kp_1 / Z_a_base;
    param.MCC_Kr_1 = MCC_Kr_1 / Z_a_base;
    param.MCC_Kp_2 = MCC_Kp_2 / Z_a_base;
    param.MCC_Kr_2 = MCC_Kr_2 / Z_a_base;
    
    param.MDV_Kp_1     = MDV_Kp_1 / 1;
    param.MDV_Ki_1     = MDV_Ki_1 / 1;
    param.MDV_C_d_1_eq = MDV_C_d_1_eq * Z_a_base;
    param.PLL_Kp       = PLL_Kp / 1;
    param.PLL_Ki       = PLL_Ki / 1;
    param.TAC_Kp_1     = TAC_Kp_1 / Z_w_base;
    param.TAC_Kr_1     = TAC_Kr_1 / Z_w_base;
    param.TAC_Kp_2     = TAC_Kp_2 / Z_w_base;
    param.TAC_Kr_2     = TAC_Kr_2 / Z_w_base;
    param.TAC_Kp_3     = TAC_Kp_3 / Z_w_base;
    param.TAC_Kr_3     = TAC_Kr_3 / Z_w_base;
    param.MPQ_Kp_1     = MPQ_Kp_1 * V_a_base;
    param.MPQ_Ki_1     = MPQ_Ki_1 * V_a_base;
    param.MPQ_Kp_2     = MPQ_Kp_2 * V_a_base;
    param.MPQ_Ki_2     = MPQ_Ki_2 * V_a_base;
    param.TPQ_Kp       = TPQ_Kp * V_w_base;
    param.TPQ_Ki       = TPQ_Ki * V_w_base;
    
    % ---------- INPUT PARAMETERS ---------- %
    
    V_g_1_ph_pk = V_g_1_base;
    I_g_1_ph_pk = I_g_1_base;
    V_g_2_ph_pk = V_g_2_base;
    I_g_2_ph_pk = I_g_2_base;
    V_a_ph_pk   = V_a_base;
    I_a_ph_pk   = I_a_base;
    V_w_ph_pk   = V_w_base;
    I_w_1_ph_pk = I_w_1_base;
    I_w_2_ph_pk = I_w_2_base;
    I_w_3_ph_pk = I_w_3_base;
    V_d_pk      = V_d_nom;
    I_d_pk      = I_d_nom;
    I_s_ph_pk   = I_a_base;
    V_d_w_pk    = V_d_w_nom;
    
    P_g_1_ref = +0.8 * S_base;
    Q_g_1_ref = -0.3 * S_base;
    P_g_2_ref = -0.8 * S_base;
    Q_g_2_ref = -0.3 * S_base;
    I_x_sum_ph_pk = 0.8 * I_g_2_ph_pk;
    
    P_f_1_ref = 0.8 * S_w_1_base;
    P_f_2_ref = 0.8 * S_w_2_base;
    P_f_3_ref = 0.8 * S_w_3_base;
    Q_f_1_ref = 0.3 * S_w_1_base;
    Q_f_2_ref = 0.3 * S_w_2_base;
    Q_f_3_ref = 0.3 * S_w_3_base;
    
    param.V_g_1_ph_pk   = V_g_1_ph_pk   / V_g_1_base;
    param.I_g_1_ph_pk   = I_g_1_ph_pk   / I_g_1_base;
    param.V_g_2_ph_pk   = V_g_2_ph_pk   / V_g_2_base;
    param.I_g_2_ph_pk   = I_g_2_ph_pk   / I_g_2_base;
    param.I_x_sum_ph_pk = I_x_sum_ph_pk / I_g_2_base;
    param.V_a_ph_pk     = V_a_ph_pk     / V_a_base;
    param.I_a_ph_pk     = I_a_ph_pk     / I_a_base;
    param.I_s_ph_pk     = I_s_ph_pk     / I_a_base;
    param.V_d_pk        = V_d_pk        / V_a_base;
    param.I_d_pk        = I_d_pk        / I_a_base;
    param.V_w_ph_pk     = V_w_ph_pk     / V_w_base;
    param.I_w_1_ph_pk   = I_w_1_ph_pk   / I_w_base;
    param.I_w_2_ph_pk   = I_w_2_ph_pk   / I_w_base;
    param.I_w_3_ph_pk   = I_w_3_ph_pk   / I_w_base;
    param.V_d_w_pk      = V_d_w_pk      / V_w_base;
    param.S_nom         = S_base        / E_base;
    param.P_g_1_ref     = P_g_1_ref     / E_base;
    param.Q_g_1_ref     = Q_g_1_ref     / E_base;
    param.P_g_2_ref     = P_g_2_ref     / E_base;
    param.Q_g_2_ref     = Q_g_2_ref     / E_base;
    param.P_f_1_ref     = P_f_1_ref     / E_base;
    param.P_f_2_ref     = P_f_2_ref     / E_base;
    param.P_f_3_ref     = P_f_3_ref     / E_base;
    param.Q_f_1_ref     = Q_f_1_ref     / E_base;
    param.Q_f_2_ref     = Q_f_2_ref     / E_base;
    param.Q_f_3_ref     = Q_f_3_ref     / E_base;
    
    % ---------- DELAY PARAMETERS ---------- %
    % setting a non-zero delay value only impacts the systems implemented
    % in such a way that they account for delays. The non-delayed systems
    % disregard the DSA.delay field.
    
    try
        T_d_MMC_us = DSA.sweep.T_d_MMC_us;
    catch
        T_d_MMC_us = 250;
        % T_d_MMC_us = 0;
    end
    T_d_TLC_us = 150;
    % T_d_TLC_us = 0;

    delay.pade.MMC = 0; % 9;
    delay.pade.TLC = 0; % 9;

    T_d_MMC = T_d_MMC_us*1e-6;
    T_d_TLC = T_d_TLC_us*1e-6; 
    
    pathname = 'C:\Users\pderua\Box Sync\research\matlab\DSAnalysis\sstms\NET\pade';
    if T_d_MMC > 0 && delay.pade.MMC > 0
        MMC_filename = [pathname '\PADE_order' num2str(delay.pade.MMC) '_' num2str(T_d_MMC_us) 'us.mat'];
        ld = load(MMC_filename);
        disp(['       MMC loading: ' MMC_filename])
        DSA.PAD.AM = ld.PADE_stsp.A;
        DSA.PAD.BM = ld.PADE_stsp.B;
        DSA.PAD.CM = ld.PADE_stsp.C;
        DSA.PAD.DM = ld.PADE_stsp.D;
    end
    if T_d_TLC > 0 && delay.pade.TLC > 0
        TLC_filename = [pathname '\PADE_order' num2str(delay.pade.TLC) '_' num2str(T_d_TLC_us) 'us.mat'];
        ld = load(TLC_filename);
        disp(['       TLC loading: ' TLC_filename])
        DSA.PAD.AT = ld.PADE_stsp.A;
        DSA.PAD.BT = ld.PADE_stsp.B;
        DSA.PAD.CT = ld.PADE_stsp.C;
        DSA.PAD.DT = ld.PADE_stsp.D;
    end
    
    delay.n_u_a_ref   = T_d_MMC;
    delay.n_u_b_ref   = T_d_MMC;
    delay.n_u_c_ref   = T_d_MMC;
    delay.n_l_a_ref   = T_d_MMC;
    delay.n_l_b_ref   = T_d_MMC;
    delay.n_l_c_ref   = T_d_MMC;
    delay.n_u_1_a_ref = T_d_MMC;
    delay.n_u_1_b_ref = T_d_MMC;
    delay.n_u_1_c_ref = T_d_MMC;
    delay.n_l_1_a_ref = T_d_MMC;
    delay.n_l_1_b_ref = T_d_MMC;
    delay.n_l_1_c_ref = T_d_MMC;
    delay.n_u_2_a_ref = T_d_MMC;
    delay.n_u_2_b_ref = T_d_MMC;
    delay.n_u_2_c_ref = T_d_MMC;
    delay.n_l_2_a_ref = T_d_MMC;
    delay.n_l_2_b_ref = T_d_MMC;
    delay.n_l_2_c_ref = T_d_MMC;
    delay.m_w_al_ref   = T_d_TLC;
    delay.m_w_be_ref   = T_d_TLC;
    delay.m_w_1_al_ref = T_d_TLC;
    delay.m_w_1_be_ref = T_d_TLC;
    delay.m_w_2_al_ref = T_d_TLC;
    delay.m_w_2_be_ref = T_d_TLC;
    delay.m_w_3_al_ref = T_d_TLC;
    delay.m_w_3_be_ref = T_d_TLC;
    
    % ---------- PARAMETERS STRUCTURES ---------- %
    
    data.f_1   = f_1;
    data.delay = delay;
    data.base  = base;
    data.param = param;
    DSA.data   = data;
end