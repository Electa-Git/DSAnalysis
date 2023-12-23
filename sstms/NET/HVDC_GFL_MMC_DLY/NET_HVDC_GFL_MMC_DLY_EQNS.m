if get_values
    omega_1      = param.omega_1;
    % common
    PLL_Kp       = param.PLL_Kp;
    PLL_Ki       = param.PLL_Ki;
    % onshore MMC
    k_MMC1       = param.k_MMC1;
    R_a_1        = param.R_a_1;
    L_a_1        = param.L_a_1;
    C_a_1        = param.C_a_1;
    R_g_1        = param.R_g_1;
    R_e_1        = param.R_e_1;
    L_e_1        = param.L_e_1;
    MAC_Kp_1     = param.MAC_Kp_1;
    MAC_Kr_1     = param.MAC_Kr_1;
    MCC_Kp_1     = param.MCC_Kp_1;
    MCC_Kr_1     = param.MCC_Kr_1;
    MPQ_Kp_1     = param.MPQ_Kp_1;
    MPQ_Ki_1     = param.MPQ_Ki_1;
    MDV_C_d_1_eq = param.MDV_C_d_1_eq;
    MDV_Kp_1     = param.MDV_Kp_1;
    MDV_Ki_1     = param.MDV_Ki_1;
    V_a_ph_pk    = param.V_a_ph_pk;
    % offshore MMC
    k_MMC2      = param.k_MMC2;
    R_a_2       = param.R_a_2;
    L_a_2       = param.L_a_2;
    C_a_2       = param.C_a_2;
    R_g_2       = param.R_g_2;
    R_e_2       = param.R_e_2;
    L_e_2       = param.L_e_2;
    MAC_Kp_2    = param.MAC_Kp_2;
    MAC_Kr_2    = param.MAC_Kr_2;
    MCC_Kp_2    = param.MCC_Kp_2;
    MCC_Kr_2    = param.MCC_Kr_2;
    MPQ_Kp_2    = param.MPQ_Kp_2;
    MPQ_Ki_2    = param.MPQ_Ki_2;
    % -- HVDC cable
    R_d         = param.R_d;
    L_d         = param.L_d;
    C_d         = param.C_d;
    G_d         = param.G_d;

    v_g_1_a       = inputs.v_g_1_a;
    v_g_1_b       = inputs.v_g_1_b;
    v_g_1_c       = inputs.v_g_1_c;
    v_g_2_a       = inputs.v_g_2_a;
    v_g_2_b       = inputs.v_g_2_b;
    v_g_2_c       = inputs.v_g_2_c;
    v_d_1_ref     = inputs.v_d_1_ref;
    q_g_1_ref     = inputs.q_g_1_ref;
    p_g_2_ref     = inputs.p_g_2_ref;
    q_g_2_ref     = inputs.q_g_2_ref;
    cos_omega_1_t = inputs.cos_omega_1_t;
    sin_omega_1_t = inputs.sin_omega_1_t;
    
    % -- onshore MMC
    i_s_1_al          = states.i_s_1_al;
    i_s_1_be          = states.i_s_1_be;
    i_c_1_a           = states.i_c_1_a;
    i_c_1_b           = states.i_c_1_b;
    i_c_1_c           = states.i_c_1_c;
    v_c_u_1_a         = states.v_c_u_1_a;
    v_c_u_1_b         = states.v_c_u_1_b;
    v_c_u_1_c         = states.v_c_u_1_c;
    v_c_l_1_a         = states.v_c_l_1_a;
    v_c_l_1_b         = states.v_c_l_1_b;
    v_c_l_1_c         = states.v_c_l_1_c;
    e_g_1_al_flt      = states.e_g_1_al_flt;
    e_g_1_be_flt      = states.e_g_1_be_flt;
    v_g_1_al_flt      = states.v_g_1_al_flt;
    v_g_1_be_flt      = states.v_g_1_be_flt;
    e_g_1_al_flt_quad = states.e_g_1_al_flt_quad;
    e_g_1_be_flt_quad = states.e_g_1_be_flt_quad;
    v_g_1_al_flt_quad = states.v_g_1_al_flt_quad;
    v_g_1_be_flt_quad = states.v_g_1_be_flt_quad;
    e_PLL_v_1         = states.e_PLL_v_1;
    theta_eps_v_1     = states.theta_eps_v_1;
    e_MV_1            = states.e_MV_1;
    e_MQ_1            = states.e_MQ_1;
    e_MAC1_1_al       = states.e_MAC1_1_al;
    e_MAC1_1_be       = states.e_MAC1_1_be;
    e_MAC2_1_al       = states.e_MAC2_1_al;
    e_MAC2_1_be       = states.e_MAC2_1_be;
    e_MCC1_1_a        = states.e_MCC1_1_a;
    e_MCC1_1_b        = states.e_MCC1_1_b;
    e_MCC1_1_c        = states.e_MCC1_1_c;
    e_MCC2_1_a        = states.e_MCC2_1_a;
    e_MCC2_1_b        = states.e_MCC2_1_b;
    e_MCC2_1_c        = states.e_MCC2_1_c;
    % -- offshore MMC
    i_s_2_al          = states.i_s_2_al;
    i_s_2_be          = states.i_s_2_be;
    i_c_2_a           = states.i_c_2_a;
    i_c_2_b           = states.i_c_2_b;
    i_c_2_c           = states.i_c_2_c;
    v_c_u_2_a         = states.v_c_u_2_a;
    v_c_u_2_b         = states.v_c_u_2_b;
    v_c_u_2_c         = states.v_c_u_2_c;
    v_c_l_2_a         = states.v_c_l_2_a;
    v_c_l_2_b         = states.v_c_l_2_b;
    v_c_l_2_c         = states.v_c_l_2_c;
    e_g_2_al_flt      = states.e_g_2_al_flt;
    e_g_2_be_flt      = states.e_g_2_be_flt;
    v_g_2_al_flt      = states.v_g_2_al_flt;
    v_g_2_be_flt      = states.v_g_2_be_flt;
    e_g_2_al_flt_quad = states.e_g_2_al_flt_quad;
    e_g_2_be_flt_quad = states.e_g_2_be_flt_quad;
    v_g_2_al_flt_quad = states.v_g_2_al_flt_quad;
    v_g_2_be_flt_quad = states.v_g_2_be_flt_quad;
    e_PLL_v_2         = states.e_PLL_v_2;
    theta_eps_v_2     = states.theta_eps_v_2;
    e_MP_2            = states.e_MP_2;
    e_MQ_2            = states.e_MQ_2;
    e_MAC1_2_al       = states.e_MAC1_2_al;
    e_MAC1_2_be       = states.e_MAC1_2_be;
    e_MAC2_2_al       = states.e_MAC2_2_al;
    e_MAC2_2_be       = states.e_MAC2_2_be;
    e_MCC1_2_a        = states.e_MCC1_2_a;
    e_MCC1_2_b        = states.e_MCC1_2_b;
    e_MCC1_2_c        = states.e_MCC1_2_c;
    e_MCC2_2_a        = states.e_MCC2_2_a;
    e_MCC2_2_b        = states.e_MCC2_2_b;
    e_MCC2_2_c        = states.e_MCC2_2_c;
    % -- HVDC cable
    v_d_u_1           = states.v_d_u_1;
    v_d_l_1           = states.v_d_l_1;
    i_z_u             = states.i_z_u;
    i_z_l             = states.i_z_l;
    v_d_u_2           = states.v_d_u_2;
    v_d_l_2           = states.v_d_l_2;
    
    n_u_1_a_ref_dlyd = delayed.n_u_1_a_ref_dlyd;
    n_u_1_b_ref_dlyd = delayed.n_u_1_b_ref_dlyd;
    n_u_1_c_ref_dlyd = delayed.n_u_1_c_ref_dlyd;
    n_l_1_a_ref_dlyd = delayed.n_l_1_a_ref_dlyd;
    n_l_1_b_ref_dlyd = delayed.n_l_1_b_ref_dlyd;
    n_l_1_c_ref_dlyd = delayed.n_l_1_c_ref_dlyd;
    n_u_2_a_ref_dlyd = delayed.n_u_2_a_ref_dlyd;
    n_u_2_b_ref_dlyd = delayed.n_u_2_b_ref_dlyd;
    n_u_2_c_ref_dlyd = delayed.n_u_2_c_ref_dlyd;
    n_l_2_a_ref_dlyd = delayed.n_l_2_a_ref_dlyd;
    n_l_2_b_ref_dlyd = delayed.n_l_2_b_ref_dlyd;
    n_l_2_c_ref_dlyd = delayed.n_l_2_c_ref_dlyd;
end

if run_algebra % carrying out pre-calculations
    % _____________________________________________________________________
    %                                                            HVDC cable
    
    v_d_1     = v_d_u_1 + v_d_l_1;
    v_d_del_1 = v_d_u_1 - v_d_l_1;
    v_d_2     = v_d_u_2 + v_d_l_2;
    v_d_del_2 = v_d_u_2 - v_d_l_2;
    
    % _____________________________________________________________________
    %                                            onshore grid-following MMC
    
    % --- Clarke transformation

    i_s_1_a = i_s_1_al;
    i_s_1_b = 1/2 * (-i_s_1_al + sqrt(3)*i_s_1_be);
    i_s_1_c = 1/2 * (-i_s_1_al - sqrt(3)*i_s_1_be);

    v_g_1_al = 2/3 * v_g_1_a - 1/3 * (v_g_1_b + v_g_1_c);
    v_g_1_be = 1/sqrt(3) * (v_g_1_b - v_g_1_c);
    
    % --- Transformer
    
    i_g_1_a = i_s_1_a / k_MMC1;
    i_g_1_b = i_s_1_b / k_MMC1;
    i_g_1_c = i_s_1_c / k_MMC1;
    
    % --- Current calculations
    
    i_u_1_a = i_c_1_a + i_s_1_a/2;
    i_u_1_b = i_c_1_b + i_s_1_b/2;
    i_u_1_c = i_c_1_c + i_s_1_c/2;
    i_l_1_a = i_c_1_a - i_s_1_a/2;
    i_l_1_b = i_c_1_b - i_s_1_b/2;
    i_l_1_c = i_c_1_c - i_s_1_c/2;
    
    i_d_u_1   = i_u_1_a + i_u_1_b + i_u_1_c;
    i_d_l_1   = i_l_1_a + i_l_1_b + i_l_1_c;
    i_d_1     = (i_d_u_1 + i_d_l_1)/2;
    i_d_del_1 = i_d_u_1 - i_d_l_1;

    % --- Positive and negative sequence extraction
    
    v_g_1_al_flt_ps = 1/2 * ( v_g_1_al_flt      - v_g_1_be_flt_quad);
    v_g_1_be_flt_ps = 1/2 * ( v_g_1_al_flt_quad + v_g_1_be_flt     );
    v_g_1_al_flt_ns = 1/2 * ( v_g_1_al_flt      + v_g_1_be_flt_quad);
    v_g_1_be_flt_ns = 1/2 * (-v_g_1_al_flt_quad + v_g_1_be_flt     );

    % --- Park transformation and PLL
    
    cos_theta_eps_1 = cos(theta_eps_v_1);
    sin_theta_eps_1 = sin(theta_eps_v_1);

    cos_theta_1 = cos_omega_1_t.*cos_theta_eps_1 - sin_omega_1_t.*sin_theta_eps_1;
    sin_theta_1 = sin_omega_1_t.*cos_theta_eps_1 + cos_omega_1_t.*sin_theta_eps_1;

    v_g_1_d_flt =   cos_theta_1.*v_g_1_al_flt_ps + sin_theta_1.*v_g_1_be_flt_ps;
    v_g_1_q_flt = - sin_theta_1.*v_g_1_al_flt_ps + cos_theta_1.*v_g_1_be_flt_ps;

    delta_theta_1 = atan2(v_g_1_q_flt, v_g_1_d_flt);
    delta_omega_1 = e_PLL_v_1 + PLL_Kp*delta_theta_1;
    
    % ---  PQ control

    p_g_1 = (v_g_1_a.*i_s_1_a + v_g_1_b.*i_s_1_b + v_g_1_c.*i_s_1_c) / k_MMC1;
    q_g_1 = 1/sqrt(3) * (v_g_1_a.*(i_s_1_c-i_s_1_b) + v_g_1_b.*(i_s_1_a-i_s_1_c) + v_g_1_c.*(i_s_1_b-i_s_1_a)) / k_MMC1;

    w_z_1_ref = MDV_C_d_1_eq/2 * v_d_1_ref.^2;
    w_z_1     = MDV_C_d_1_eq/2 * v_d_1.^2;
    
    i_s_1_d_ref = (e_MV_1 - MDV_Kp_1 * (w_z_1_ref - w_z_1)) ./ (3/2 * V_a_ph_pk);    
    i_s_1_q_ref = e_MQ_1 - MPQ_Kp_1 * (q_g_1_ref - q_g_1);

    i_s_1_al_ref = cos_theta_1.*i_s_1_d_ref - sin_theta_1.*i_s_1_q_ref;
    i_s_1_be_ref = sin_theta_1.*i_s_1_d_ref + cos_theta_1.*i_s_1_q_ref;

    % --- AC control

    v_s_1_al_ref = v_g_1_al_flt/k_MMC1 + e_MAC2_1_al + MAC_Kp_1 * (i_s_1_al_ref - i_s_1_al);
    v_s_1_be_ref = v_g_1_be_flt/k_MMC1 + e_MAC2_1_be + MAC_Kp_1 * (i_s_1_be_ref - i_s_1_be);

    v_s_1_a_ref = v_s_1_al_ref;
    v_s_1_b_ref = 1/2 * (-v_s_1_al_ref + sqrt(3)*v_s_1_be_ref);
    v_s_1_c_ref = 1/2 * (-v_s_1_al_ref - sqrt(3)*v_s_1_be_ref);
    
    % ---  CC control

    p_1_loss = (i_s_1_a.*i_s_1_a + i_s_1_b.*i_s_1_b + i_s_1_c.*i_s_1_c) * R_g_1;
    i_c_1_ref = (p_g_1 + p_1_loss) ./ (3*v_d_1);

    v_c_1_a_ref = v_d_1/2 + e_MCC2_1_a - MCC_Kp_1 * (i_c_1_ref - i_c_1_a);
    v_c_1_b_ref = v_d_1/2 + e_MCC2_1_b - MCC_Kp_1 * (i_c_1_ref - i_c_1_b);
    v_c_1_c_ref = v_d_1/2 + e_MCC2_1_c - MCC_Kp_1 * (i_c_1_ref - i_c_1_c);

    % --- modulation indices

    n_u_1_a_ref = (-v_s_1_a_ref + v_c_1_a_ref) ./ v_d_1;
    n_u_1_b_ref = (-v_s_1_b_ref + v_c_1_b_ref) ./ v_d_1;
    n_u_1_c_ref = (-v_s_1_c_ref + v_c_1_c_ref) ./ v_d_1;
    n_l_1_a_ref = ( v_s_1_a_ref + v_c_1_a_ref) ./ v_d_1;
    n_l_1_b_ref = ( v_s_1_b_ref + v_c_1_b_ref) ./ v_d_1;
    n_l_1_c_ref = ( v_s_1_c_ref + v_c_1_c_ref) ./ v_d_1;
    
    if is_delayed
        n_u_1_a = n_u_1_a_ref_dlyd;
        n_u_1_b = n_u_1_b_ref_dlyd;
        n_u_1_c = n_u_1_c_ref_dlyd;
        n_l_1_a = n_l_1_a_ref_dlyd;
        n_l_1_b = n_l_1_b_ref_dlyd;
        n_l_1_c = n_l_1_c_ref_dlyd;
    else
        n_u_1_a = n_u_1_a_ref;
        n_u_1_b = n_u_1_b_ref;
        n_u_1_c = n_u_1_c_ref;
        n_l_1_a = n_l_1_a_ref;
        n_l_1_b = n_l_1_b_ref;
        n_l_1_c = n_l_1_c_ref;
    end
    
    % --- inserted and driving voltages
    
    v_u_1_a = n_u_1_a.*v_c_u_1_a;
    v_u_1_b = n_u_1_b.*v_c_u_1_b;
    v_u_1_c = n_u_1_c.*v_c_u_1_c;
    v_l_1_a = n_l_1_a.*v_c_l_1_a;
    v_l_1_b = n_l_1_b.*v_c_l_1_b;
    v_l_1_c = n_l_1_c.*v_c_l_1_c;
    
    v_s_1_a = 1/2 * (-v_u_1_a + v_l_1_a);
    v_s_1_b = 1/2 * (-v_u_1_b + v_l_1_b);
    v_s_1_c = 1/2 * (-v_u_1_c + v_l_1_c);
    v_c_1_a = 1/2 * ( v_u_1_a + v_l_1_a);
    v_c_1_b = 1/2 * ( v_u_1_b + v_l_1_b);
    v_c_1_c = 1/2 * ( v_u_1_c + v_l_1_c);
    
    v_s_1_al = 2/3 * v_s_1_a - 1/3 * (v_s_1_b + v_s_1_c);
    v_s_1_be = 1/sqrt(3) * (v_s_1_b - v_s_1_c);

    % _____________________________________________________________________
    %                                           offshore grid-following MMC

    % --- Clarke transformations
    
    i_s_2_a = i_s_2_al;
    i_s_2_b = 1/2 * (-i_s_2_al + sqrt(3)*i_s_2_be);
    i_s_2_c = 1/2 * (-i_s_2_al - sqrt(3)*i_s_2_be);

    v_g_2_al = 2/3 * v_g_2_a - 1/3 * (v_g_2_b + v_g_2_c);
    v_g_2_be = 1/sqrt(3) * (v_g_2_b - v_g_2_c);
    
    % --- Transformer
    
    i_g_2_a = i_s_2_a / k_MMC2;
    i_g_2_b = i_s_2_b / k_MMC2;
    i_g_2_c = i_s_2_c / k_MMC2;

    % --- Current calculations
    
    i_u_2_a = i_c_2_a + i_s_2_a/2;
    i_u_2_b = i_c_2_b + i_s_2_b/2;
    i_u_2_c = i_c_2_c + i_s_2_c/2;
    i_l_2_a = i_c_2_a - i_s_2_a/2;
    i_l_2_b = i_c_2_b - i_s_2_b/2;
    i_l_2_c = i_c_2_c - i_s_2_c/2;

    i_d_u_2   = i_u_2_a + i_u_2_b + i_u_2_c;
    i_d_l_2   = i_l_2_a + i_l_2_b + i_l_2_c;
    i_d_2     = (i_d_u_2 + i_d_l_2)/2;
    i_d_del_2 = i_d_u_2 - i_d_l_2;
    
    % --- Positive and negative sequence extraction
    
    v_g_2_al_flt_ps = 1/2 * ( v_g_2_al_flt      - v_g_2_be_flt_quad);
    v_g_2_be_flt_ps = 1/2 * ( v_g_2_al_flt_quad + v_g_2_be_flt     );
    v_g_2_al_flt_ns = 1/2 * ( v_g_2_al_flt      + v_g_2_be_flt_quad);
    v_g_2_be_flt_ns = 1/2 * (-v_g_2_al_flt_quad + v_g_2_be_flt     );

    v_g_2_a_flt = v_g_2_al_flt;
    v_g_2_b_flt = 1/2 * (-v_g_2_al_flt + sqrt(3)*v_g_2_be_flt);
    v_g_2_c_flt = 1/2 * (-v_g_2_al_flt - sqrt(3)*v_g_2_be_flt);
    
    % --- Park transformation and PLL
    
    cos_theta_eps_2 = cos(theta_eps_v_2);
    sin_theta_eps_2 = sin(theta_eps_v_2);

    cos_theta_2 = cos_omega_1_t.*cos_theta_eps_2 - sin_omega_1_t.*sin_theta_eps_2;
    sin_theta_2 = sin_omega_1_t.*cos_theta_eps_2 + cos_omega_1_t.*sin_theta_eps_2;

    v_g_2_d_flt =   cos_theta_2.*v_g_2_al_flt_ps + sin_theta_2.*v_g_2_be_flt_ps;
    v_g_2_q_flt = - sin_theta_2.*v_g_2_al_flt_ps + cos_theta_2.*v_g_2_be_flt_ps;

    delta_theta_2 = atan2(v_g_2_q_flt, v_g_2_d_flt);
    delta_omega_2 = e_PLL_v_2 + PLL_Kp*delta_theta_2;
    
    % ---  PQ control

    p_g_2 = (v_g_2_a.*i_s_2_a + v_g_2_b.*i_s_2_b + v_g_2_c.*i_s_2_c) / k_MMC2;
    q_g_2 = 1/sqrt(3) * (v_g_2_a.*(i_s_2_c-i_s_2_b) + v_g_2_b.*(i_s_2_a-i_s_2_c) + v_g_2_c.*(i_s_2_b-i_s_2_a)) / k_MMC2;
    
    i_s_2_d_ref = e_MP_2 + MPQ_Kp_2 * (p_g_2_ref - p_g_2);
    i_s_2_q_ref = e_MQ_2 - MPQ_Kp_2 * (q_g_2_ref - q_g_2);
    
    i_s_2_al_ref = cos_theta_2.*i_s_2_d_ref - sin_theta_2.*i_s_2_q_ref;
    i_s_2_be_ref = sin_theta_2.*i_s_2_d_ref + cos_theta_2.*i_s_2_q_ref;

    i_s_2_a_ref = i_s_2_al_ref;
    i_s_2_b_ref = 1/2 * (-i_s_2_al_ref + sqrt(3)*i_s_2_be_ref);
    i_s_2_c_ref = 1/2 * (-i_s_2_al_ref - sqrt(3)*i_s_2_be_ref);
    
    % --- alternating current control
    
    v_s_2_al_ref = v_g_2_al_flt + e_MAC2_2_al + MAC_Kp_2 * (i_s_2_al_ref - i_s_2_al);
    v_s_2_be_ref = v_g_2_be_flt + e_MAC2_2_be + MAC_Kp_2 * (i_s_2_be_ref - i_s_2_be);

    v_s_2_a_ref = v_s_2_al_ref;
    v_s_2_b_ref = 1/2 * (-v_s_2_al_ref + sqrt(3)*v_s_2_be_ref);
    v_s_2_c_ref = 1/2 * (-v_s_2_al_ref - sqrt(3)*v_s_2_be_ref);
    
    % ---  circulating current control

    p_2_loss = (i_s_2_a.*i_s_2_a + i_s_2_b.*i_s_2_b + i_s_2_c.*i_s_2_c) * R_g_2;
    i_c_2_ref = (p_g_2 + p_2_loss) ./ (3*v_d_2);
    
    v_c_2_a_ref = v_d_2/2 + e_MCC2_2_a - MCC_Kp_2 * (i_c_2_ref - i_c_2_a);
    v_c_2_b_ref = v_d_2/2 + e_MCC2_2_b - MCC_Kp_2 * (i_c_2_ref - i_c_2_b);
    v_c_2_c_ref = v_d_2/2 + e_MCC2_2_c - MCC_Kp_2 * (i_c_2_ref - i_c_2_c);

    % --- modulation indices

    n_u_2_a_ref = (-v_s_2_a_ref + v_c_2_a_ref) ./ v_d_2;
    n_u_2_b_ref = (-v_s_2_b_ref + v_c_2_b_ref) ./ v_d_2;
    n_u_2_c_ref = (-v_s_2_c_ref + v_c_2_c_ref) ./ v_d_2;
    n_l_2_a_ref = ( v_s_2_a_ref + v_c_2_a_ref) ./ v_d_2;
    n_l_2_b_ref = ( v_s_2_b_ref + v_c_2_b_ref) ./ v_d_2;
    n_l_2_c_ref = ( v_s_2_c_ref + v_c_2_c_ref) ./ v_d_2;
    
    if is_delayed
        n_u_2_a = n_u_2_a_ref_dlyd;
        n_u_2_b = n_u_2_b_ref_dlyd;
        n_u_2_c = n_u_2_c_ref_dlyd;
        n_l_2_a = n_l_2_a_ref_dlyd;
        n_l_2_b = n_l_2_b_ref_dlyd;
        n_l_2_c = n_l_2_c_ref_dlyd;
    else
        n_u_2_a = n_u_2_a_ref;
        n_u_2_b = n_u_2_b_ref;
        n_u_2_c = n_u_2_c_ref;
        n_l_2_a = n_l_2_a_ref;
        n_l_2_b = n_l_2_b_ref;
        n_l_2_c = n_l_2_c_ref;
    end
    
    % --- inserted and driving voltages
    
    v_u_2_a = n_u_2_a.*v_c_u_2_a;
    v_u_2_b = n_u_2_b.*v_c_u_2_b;
    v_u_2_c = n_u_2_c.*v_c_u_2_c;
    v_l_2_a = n_l_2_a.*v_c_l_2_a;
    v_l_2_b = n_l_2_b.*v_c_l_2_b;
    v_l_2_c = n_l_2_c.*v_c_l_2_c;
    
    v_s_2_a = 1/2 * (-v_u_2_a + v_l_2_a);
    v_s_2_b = 1/2 * (-v_u_2_b + v_l_2_b);
    v_s_2_c = 1/2 * (-v_u_2_c + v_l_2_c);
    v_c_2_a = 1/2 * ( v_u_2_a + v_l_2_a);
    v_c_2_b = 1/2 * ( v_u_2_b + v_l_2_b);
    v_c_2_c = 1/2 * ( v_u_2_c + v_l_2_c);
    
    v_s_2_al = 2/3 * v_s_2_a - 1/3 * (v_s_2_b + v_s_2_c);
    v_s_2_be = 1/sqrt(3) * (v_s_2_b - v_s_2_c);
end

if run_dxdt % defining dxdt = f(x,u)
    dxdt = [
        % _____________________________________________________________________
        %                                            onshore grid-following MMC
        1/L_e_1 * (-R_e_1*i_s_1_al + v_s_1_al - v_g_1_al/k_MMC1);        % i_s_1_al
        1/L_e_1 * (-R_e_1*i_s_1_be + v_s_1_be - v_g_1_be/k_MMC1);        % i_s_1_be
        1/L_a_1 * (-R_a_1*i_c_1_a - v_c_1_a + v_d_1/2);                  % i_c_1_a
        1/L_a_1 * (-R_a_1*i_c_1_b - v_c_1_b + v_d_1/2);                  % i_c_1_b
        1/L_a_1 * (-R_a_1*i_c_1_c - v_c_1_c + v_d_1/2);                  % i_c_1_c
        1/C_a_1 * n_u_1_a.*i_u_1_a;                                      % v_c_u_1_a
        1/C_a_1 * n_u_1_b.*i_u_1_b;                                      % v_c_u_1_b
        1/C_a_1 * n_u_1_c.*i_u_1_c;                                      % v_c_u_1_c
        1/C_a_1 * n_l_1_a.*i_l_1_a;                                      % v_c_l_1_a
        1/C_a_1 * n_l_1_b.*i_l_1_b;                                      % v_c_l_1_b
        1/C_a_1 * n_l_1_c.*i_l_1_c;                                      % v_c_l_1_c
        
        -omega_1 * v_g_1_al_flt;                                         % e_g_1_al_flt
        -omega_1 * v_g_1_be_flt;                                         % e_g_1_be_flt
        +omega_1 * (e_g_1_al_flt + v_g_1_al - v_g_1_al_flt);             % v_g_1_al_flt
        +omega_1 * (e_g_1_be_flt + v_g_1_be - v_g_1_be_flt);             % v_g_1_be_flt
        +omega_1 * (v_g_1_al - v_g_1_al_flt_quad);                       % e_g_1_al_flt_quad
        +omega_1 * (v_g_1_be - v_g_1_be_flt_quad);                       % e_g_1_be_flt_quad
        +omega_1 * (e_g_1_al_flt_quad - v_g_1_al_flt_quad);              % v_g_1_al_flt_quad
        +omega_1 * (e_g_1_be_flt_quad - v_g_1_be_flt_quad);              % v_g_1_be_flt_quad
        
        PLL_Ki * delta_theta_1;                                          % e_PLL_v_1
        delta_omega_1;                                                   % theta_eps_v_1
        -MDV_Ki_1 * (w_z_1_ref - w_z_1);                                 % e_MV_1
        -MPQ_Ki_1 * (q_g_1_ref - q_g_1);                                 % e_MQ_1
        -omega_1 * e_MAC2_1_al;                                          % e_MAC1_1_al
        -omega_1 * e_MAC2_1_be;                                          % e_MAC1_1_be
        +omega_1 * e_MAC1_1_al + MAC_Kr_1 * (i_s_1_al_ref - i_s_1_al);   % e_MAC2_1_al
        +omega_1 * e_MAC1_1_be + MAC_Kr_1 * (i_s_1_be_ref - i_s_1_be);   % e_MAC2_1_be
        -2*omega_1 * e_MCC2_1_a;                                         % e_MCC1_1_a
        -2*omega_1 * e_MCC2_1_b;                                         % e_MCC1_1_b
        -2*omega_1 * e_MCC2_1_c;                                         % e_MCC1_1_c
        +2*omega_1 * e_MCC1_1_a - MCC_Kr_1 * (i_c_1_ref - i_c_1_a);      % e_MCC2_1_a
        +2*omega_1 * e_MCC1_1_b - MCC_Kr_1 * (i_c_1_ref - i_c_1_b);      % e_MCC2_1_b
        +2*omega_1 * e_MCC1_1_c - MCC_Kr_1 * (i_c_1_ref - i_c_1_c);      % e_MCC2_1_c
        
        % _____________________________________________________________________
        %                                           offshore grid-following MMC
        1/L_e_2 * (-R_e_2*i_s_2_al + v_s_2_al - v_g_2_al/k_MMC2);        % i_s_2_al
        1/L_e_2 * (-R_e_2*i_s_2_be + v_s_2_be - v_g_2_be/k_MMC2);        % i_s_2_be
        1/L_a_2 * (-R_a_2*i_c_2_a - v_c_2_a + v_d_2/2);                  % i_c_2_a
        1/L_a_2 * (-R_a_2*i_c_2_b - v_c_2_b + v_d_2/2);                  % i_c_2_b
        1/L_a_2 * (-R_a_2*i_c_2_c - v_c_2_c + v_d_2/2);                  % i_c_2_c
        1/C_a_2 * n_u_2_a.*i_u_2_a;                                      % v_c_u_2_a
        1/C_a_2 * n_u_2_b.*i_u_2_b;                                      % v_c_u_2_b
        1/C_a_2 * n_u_2_c.*i_u_2_c;                                      % v_c_u_2_c
        1/C_a_2 * n_l_2_a.*i_l_2_a;                                      % v_c_l_2_a
        1/C_a_2 * n_l_2_b.*i_l_2_b;                                      % v_c_l_2_b
        1/C_a_2 * n_l_2_c.*i_l_2_c;                                      % v_c_l_2_c
        
        -omega_1 * v_g_2_al_flt;                                         % e_g_2_al_flt
        -omega_1 * v_g_2_be_flt;                                         % e_g_2_be_flt
        +omega_1 * (e_g_2_al_flt + v_g_2_al - v_g_2_al_flt);             % v_g_2_al_flt
        +omega_1 * (e_g_2_be_flt + v_g_2_be - v_g_2_be_flt);             % v_g_2_be_flt
        +omega_1 * (v_g_2_al - v_g_2_al_flt_quad);                       % e_g_2_al_flt_quad
        +omega_1 * (v_g_2_be - v_g_2_be_flt_quad);                       % e_g_2_be_flt_quad
        +omega_1 * (e_g_2_al_flt_quad - v_g_2_al_flt_quad);              % v_g_2_al_flt_quad
        +omega_1 * (e_g_2_be_flt_quad - v_g_2_be_flt_quad);              % v_g_2_be_flt_quad
        
        PLL_Ki * delta_theta_2;                                          % e_PLL_v_2
        delta_omega_2;                                                   % theta_eps_v_2
        + MPQ_Ki_2 * (p_g_2_ref - p_g_2);                                % e_MP_2
        - MPQ_Ki_2 * (q_g_2_ref - q_g_2);                                % e_MQ_2
        -omega_1 * e_MAC2_2_al;                                          % e_MAC1_2_al
        -omega_1 * e_MAC2_2_be;                                          % e_MAC1_2_be
        +omega_1 * e_MAC1_2_al + MAC_Kr_2 * (i_s_2_al_ref - i_s_2_al);   % e_MAC2_2_al
        +omega_1 * e_MAC1_2_be + MAC_Kr_2 * (i_s_2_be_ref - i_s_2_be);   % e_MAC2_2_be
        -2*omega_1 * e_MCC2_2_a;                                         % e_CC1_2_a
        -2*omega_1 * e_MCC2_2_b;                                         % e_CC1_2_b
        -2*omega_1 * e_MCC2_2_c;                                         % e_CC1_2_c
        +2*omega_1 * e_MCC1_2_a - MCC_Kr_2 * (i_c_2_ref - i_c_2_a);      % e_CC2_2_a
        +2*omega_1 * e_MCC1_2_b - MCC_Kr_2 * (i_c_2_ref - i_c_2_b);      % e_CC2_2_b
        +2*omega_1 * e_MCC1_2_c - MCC_Kr_2 * (i_c_2_ref - i_c_2_c);      % e_CC2_2_c
        
        % _____________________________________________________________________
        %                                                            HVDC cable
        1/(C_d/2) * (-G_d*v_d_u_1 + i_z_u - i_d_u_1);                  % v_d_u_1
        1/(C_d/2) * (-G_d*v_d_l_1 + i_z_l - i_d_l_1);                  % v_d_l_1
        1/L_d     * (-R_d*i_z_u + v_d_u_2 - v_d_u_1)                   % i_z_u
        1/L_d     * (-R_d*i_z_l + v_d_l_2 - v_d_l_1)                   % i_z_l
        1/(C_d/2) * (-G_d*v_d_u_2 - i_z_u - i_d_u_2);                  % v_d_u_2
        1/(C_d/2) * (-G_d*v_d_l_2 - i_z_l - i_d_l_2);                  % v_d_l_2
        ];
end