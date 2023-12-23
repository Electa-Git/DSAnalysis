if get_values
    omega_1  = param.omega_1;
    R_f_1    = param.R_f_1;
    L_f_1    = param.L_f_1;
    C_f_1    = param.C_f_1;
    R_t_1    = param.R_t_1;
    L_t_1    = param.L_t_1;
    R_x_1    = param.R_x_1;
    L_x_1    = param.L_x_1;
    C_x_1    = param.C_x_1;
    R_f_2    = param.R_f_2;
    L_f_2    = param.L_f_2;
    C_f_2    = param.C_f_2;
    R_t_2    = param.R_t_2;
    L_t_2    = param.L_t_2;
    R_x_2    = param.R_x_2;
    L_x_2    = param.L_x_2;
    C_x_2    = param.C_x_2;
    R_f_3    = param.R_f_3;
    L_f_3    = param.L_f_3;
    C_f_3    = param.C_f_3;
    R_t_3    = param.R_t_3;
    L_t_3    = param.L_t_3;
    R_x_3    = param.R_x_3;
    L_x_3    = param.L_x_3;
    C_x_3    = param.C_x_3;
    k_TLC    = param.k_TLC;
    PLL_Kp   = param.PLL_Kp;
    PLL_Ki   = param.PLL_Ki;
    TAC_Kp_1 = param.TAC_Kp_1;
    TAC_Kr_1 = param.TAC_Kr_1;
    TAC_Kp_2 = param.TAC_Kp_2;
    TAC_Kr_2 = param.TAC_Kr_2;
    TAC_Kp_3 = param.TAC_Kp_3;
    TAC_Kr_3 = param.TAC_Kr_3;
    TPQ_Kp   = param.TPQ_Kp;
    TPQ_Ki   = param.TPQ_Ki;

    v_g_2_a       = inputs.v_g_2_a;
    v_g_2_b       = inputs.v_g_2_b;
    v_g_2_c       = inputs.v_g_2_c;
    p_f_1_ref     = inputs.p_f_1_ref;
    q_f_1_ref     = inputs.q_f_1_ref;
    p_f_2_ref     = inputs.p_f_2_ref;
    q_f_2_ref     = inputs.q_f_2_ref;
    p_f_3_ref     = inputs.p_f_3_ref;
    q_f_3_ref     = inputs.q_f_3_ref;
    v_d_w_1       = inputs.v_d_w_1;
    v_d_w_2       = inputs.v_d_w_2;
    v_d_w_3       = inputs.v_d_w_3;
    cos_omega_1_t = inputs.cos_omega_1_t;
    sin_omega_1_t = inputs.sin_omega_1_t;
    
    % TCL 1 and collection cable 1
    i_f_1_al          = states.i_f_1_al         ;
    i_f_1_be          = states.i_f_1_be         ;
    v_f_1_al          = states.v_f_1_al         ;
    v_f_1_be          = states.v_f_1_be         ;
    i_t_1_al          = states.i_t_1_al         ;
    i_t_1_be          = states.i_t_1_be         ;
    v_x_1_al          = states.v_x_1_al         ;
    v_x_1_be          = states.v_x_1_be         ;
    i_x_1_al          = states.i_x_1_al         ;
    i_x_1_be          = states.i_x_1_be         ;
    e_f_1_al_flt      = states.e_f_1_al_flt     ;
    e_f_1_be_flt      = states.e_f_1_be_flt     ;
    v_f_1_al_flt      = states.v_f_1_al_flt     ;
    v_f_1_be_flt      = states.v_f_1_be_flt     ;
    e_f_1_al_flt_quad = states.e_f_1_al_flt_quad;
    e_f_1_be_flt_quad = states.e_f_1_be_flt_quad;
    v_f_1_al_flt_quad = states.v_f_1_al_flt_quad;
    v_f_1_be_flt_quad = states.v_f_1_be_flt_quad;
    e_PLL_w_1         = states.e_PLL_w_1        ;
    theta_eps_w_1     = states.theta_eps_w_1    ;
    e_TPQ_1_d         = states.e_TPQ_1_d        ;
    e_TPQ_1_q         = states.e_TPQ_1_q        ;
    e_TAC1_1_al       = states.e_TAC1_1_al      ;
    e_TAC1_1_be       = states.e_TAC1_1_be      ;
    e_TAC2_1_al       = states.e_TAC2_1_al      ;
    e_TAC2_1_be       = states.e_TAC2_1_be      ;
    % TCL 2 and collection cable 2
    i_f_2_al          = states.i_f_2_al         ;
    i_f_2_be          = states.i_f_2_be         ;
    v_f_2_al          = states.v_f_2_al         ;
    v_f_2_be          = states.v_f_2_be         ;
    i_t_2_al          = states.i_t_2_al         ;
    i_t_2_be          = states.i_t_2_be         ;
    v_x_2_al          = states.v_x_2_al         ;
    v_x_2_be          = states.v_x_2_be         ;
    i_x_2_al          = states.i_x_2_al         ;
    i_x_2_be          = states.i_x_2_be         ;
    e_f_2_al_flt      = states.e_f_2_al_flt     ;
    e_f_2_be_flt      = states.e_f_2_be_flt     ;
    v_f_2_al_flt      = states.v_f_2_al_flt     ;
    v_f_2_be_flt      = states.v_f_2_be_flt     ;
    e_f_2_al_flt_quad = states.e_f_2_al_flt_quad;
    e_f_2_be_flt_quad = states.e_f_2_be_flt_quad;
    v_f_2_al_flt_quad = states.v_f_2_al_flt_quad;
    v_f_2_be_flt_quad = states.v_f_2_be_flt_quad;
    e_PLL_w_2         = states.e_PLL_w_2        ;
    theta_eps_w_2     = states.theta_eps_w_2    ;
    e_TPQ_2_d         = states.e_TPQ_2_d        ;
    e_TPQ_2_q         = states.e_TPQ_2_q        ;
    e_TAC1_2_al       = states.e_TAC1_2_al      ;
    e_TAC1_2_be       = states.e_TAC1_2_be      ;
    e_TAC2_2_al       = states.e_TAC2_2_al      ;
    e_TAC2_2_be       = states.e_TAC2_2_be      ;
    % TCL 3 and collection cable 3
    i_f_3_al          = states.i_f_3_al         ;
    i_f_3_be          = states.i_f_3_be         ;
    v_f_3_al          = states.v_f_3_al         ;
    v_f_3_be          = states.v_f_3_be         ;
    i_t_3_al          = states.i_t_3_al         ;
    i_t_3_be          = states.i_t_3_be         ;
    v_x_3_al          = states.v_x_3_al         ;
    v_x_3_be          = states.v_x_3_be         ;
    i_x_3_al          = states.i_x_3_al         ;
    i_x_3_be          = states.i_x_3_be         ;
    e_f_3_al_flt      = states.e_f_3_al_flt     ;
    e_f_3_be_flt      = states.e_f_3_be_flt     ;
    v_f_3_al_flt      = states.v_f_3_al_flt     ;
    v_f_3_be_flt      = states.v_f_3_be_flt     ;
    e_f_3_al_flt_quad = states.e_f_3_al_flt_quad;
    e_f_3_be_flt_quad = states.e_f_3_be_flt_quad;
    v_f_3_al_flt_quad = states.v_f_3_al_flt_quad;
    v_f_3_be_flt_quad = states.v_f_3_be_flt_quad;
    e_PLL_w_3         = states.e_PLL_w_3        ;
    theta_eps_w_3     = states.theta_eps_w_3    ;
    e_TPQ_3_d         = states.e_TPQ_3_d        ;
    e_TPQ_3_q         = states.e_TPQ_3_q        ;
    e_TAC1_3_al       = states.e_TAC1_3_al      ;
    e_TAC1_3_be       = states.e_TAC1_3_be      ;
    e_TAC2_3_al       = states.e_TAC2_3_al      ;
    e_TAC2_3_be       = states.e_TAC2_3_be      ;
    
    m_w_1_al_ref_dlyd = delayed.m_w_1_al_ref_dlyd;
    m_w_1_be_ref_dlyd = delayed.m_w_1_be_ref_dlyd;
    m_w_2_al_ref_dlyd = delayed.m_w_2_al_ref_dlyd;
    m_w_2_be_ref_dlyd = delayed.m_w_2_be_ref_dlyd;
    m_w_3_al_ref_dlyd = delayed.m_w_3_al_ref_dlyd;
    m_w_3_be_ref_dlyd = delayed.m_w_3_be_ref_dlyd;
end

if run_algebra
    
    % _____________________________________________________________________
    %                                                                 TLC 1
    % --- Clarke transformations
    
    i_f_1_a = i_f_1_al;
    i_f_1_b = 1/2 * (-i_f_1_al + sqrt(3)*i_f_1_be);
    i_f_1_c = 1/2 * (-i_f_1_al - sqrt(3)*i_f_1_be);
    
    i_x_1_a = i_x_1_al;
    i_x_1_b = 1/2 * (-i_x_1_al + sqrt(3)*i_x_1_be);
    i_x_1_c = 1/2 * (-i_x_1_al - sqrt(3)*i_x_1_be);
    
    v_f_1_a = v_f_1_al;
    v_f_1_b = 1/2 * (-v_f_1_al + sqrt(3)*v_f_1_be);
    v_f_1_c = 1/2 * (-v_f_1_al - sqrt(3)*v_f_1_be);
    
    % --- Transformer

    v_t_1_al = v_x_1_al / k_TLC;
    v_t_1_be = v_x_1_be / k_TLC;
    
    i_y_1_al = i_t_1_al / k_TLC;
    i_y_1_be = i_t_1_be / k_TLC;

    % --- Positive and negative sequence extraction
    
    v_f_1_al_flt_ps = 1/2 * ( v_f_1_al_flt      - v_f_1_be_flt_quad);
    v_f_1_be_flt_ps = 1/2 * ( v_f_1_al_flt_quad + v_f_1_be_flt     );
    v_f_1_al_flt_ns = 1/2 * ( v_f_1_al_flt      + v_f_1_be_flt_quad);
    v_f_1_be_flt_ns = 1/2 * (-v_f_1_al_flt_quad + v_f_1_be_flt     );
    
    % --- Park transformation and PLL
    
    cos_theta_eps_w_1 = cos(theta_eps_w_1);
    sin_theta_eps_w_1 = sin(theta_eps_w_1);

    cos_theta_w_1 = cos_omega_1_t.*cos_theta_eps_w_1 - sin_omega_1_t.*sin_theta_eps_w_1;
    sin_theta_w_1 = sin_omega_1_t.*cos_theta_eps_w_1 + cos_omega_1_t.*sin_theta_eps_w_1;

    % the PLL synchronises on the filtered positive-sequence voltage:
    v_f_1_d_flt_ps =   cos_theta_w_1.*v_f_1_al_flt_ps + sin_theta_w_1.*v_f_1_be_flt_ps;
    v_f_1_q_flt_ps = - sin_theta_w_1.*v_f_1_al_flt_ps + cos_theta_w_1.*v_f_1_be_flt_ps;
    
    delta_theta_w_1 = atan2(v_f_1_q_flt_ps, v_f_1_d_flt_ps);
    delta_omega_w_1 = e_PLL_w_1 + PLL_Kp*delta_theta_w_1;

    % --- PQ control
    
    p_f_1 = v_f_1_a.*i_f_1_a + v_f_1_b.*i_f_1_b + v_f_1_c.*i_f_1_c;
    q_f_1 = 1/sqrt(3) * (v_f_1_a.*(i_f_1_c-i_f_1_b) + v_f_1_b.*(i_f_1_a-i_f_1_c) + v_f_1_c.*(i_f_1_b-i_f_1_a));

    i_f_1_d_ref = e_TPQ_1_d + TPQ_Kp * (p_f_1_ref - p_f_1);
    i_f_1_q_ref = e_TPQ_1_q - TPQ_Kp * (q_f_1_ref - q_f_1);

    i_f_1_al_ref = cos_theta_w_1.*i_f_1_d_ref - sin_theta_w_1.*i_f_1_q_ref;
    i_f_1_be_ref = sin_theta_w_1.*i_f_1_d_ref + cos_theta_w_1.*i_f_1_q_ref;

    % --- AC control

    v_w_1_al_ref = v_f_1_al_flt + e_TAC2_1_al + TAC_Kp_1 * (i_f_1_al_ref - i_f_1_al);
    v_w_1_be_ref = v_f_1_be_flt + e_TAC2_1_be + TAC_Kp_1 * (i_f_1_be_ref - i_f_1_be);

    % --- Modulation indices
    
    m_w_1_al_ref = v_w_1_al_ref ./ (v_d_w_1/2);
    m_w_1_be_ref = v_w_1_be_ref ./ (v_d_w_1/2);

    if is_delayed
        m_w_1_al = m_w_1_al_ref_dlyd;
        m_w_1_be = m_w_1_be_ref_dlyd;
    else
        m_w_1_al = m_w_1_al_ref;
        m_w_1_be = m_w_1_be_ref;
    end

    % --- voltage modulation

    v_w_1_al = v_d_w_1/2 .* m_w_1_al;
    v_w_1_be = v_d_w_1/2 .* m_w_1_be;
    
    % --- direct current
    i_d_w_1 = 3/4 * (m_w_1_al.*i_f_1_al + m_w_1_be.*i_f_1_be);

    % _____________________________________________________________________
    %                                                                 TLC_2
    % --- Clarke transformations
    
    i_f_2_a = i_f_2_al;
    i_f_2_b = 1/2 * (-i_f_2_al + sqrt(3)*i_f_2_be);
    i_f_2_c = 1/2 * (-i_f_2_al - sqrt(3)*i_f_2_be);
    
    i_x_2_a = i_x_2_al;
    i_x_2_b = 1/2 * (-i_x_2_al + sqrt(3)*i_x_2_be);
    i_x_2_c = 1/2 * (-i_x_2_al - sqrt(3)*i_x_2_be);
    
    v_f_2_a = v_f_2_al;
    v_f_2_b = 1/2 * (-v_f_2_al + sqrt(3)*v_f_2_be);
    v_f_2_c = 1/2 * (-v_f_2_al - sqrt(3)*v_f_2_be);
    
    % --- Transformer

    v_t_2_al = v_x_2_al / k_TLC;
    v_t_2_be = v_x_2_be / k_TLC;
    
    i_y_2_al = i_t_2_al / k_TLC;
    i_y_2_be = i_t_2_be / k_TLC;

    % --- Positive and negative sequence extraction
    
    v_f_2_al_flt_ps = 1/2 * ( v_f_2_al_flt      - v_f_2_be_flt_quad);
    v_f_2_be_flt_ps = 1/2 * ( v_f_2_al_flt_quad + v_f_2_be_flt     );
    v_f_2_al_flt_ns = 1/2 * ( v_f_2_al_flt      + v_f_2_be_flt_quad);
    v_f_2_be_flt_ns = 1/2 * (-v_f_2_al_flt_quad + v_f_2_be_flt     );
    
    % --- Park transformation and PLL
    
    cos_theta_eps_w_2 = cos(theta_eps_w_2);
    sin_theta_eps_w_2 = sin(theta_eps_w_2);

    cos_theta_w_2 = cos_omega_1_t.*cos_theta_eps_w_2 - sin_omega_1_t.*sin_theta_eps_w_2;
    sin_theta_w_2 = sin_omega_1_t.*cos_theta_eps_w_2 + cos_omega_1_t.*sin_theta_eps_w_2;

    % the PLL synchronises on the filtered positive-sequence voltage:
    v_f_2_d_flt_ps =   cos_theta_w_2.*v_f_2_al_flt_ps + sin_theta_w_2.*v_f_2_be_flt_ps;
    v_f_2_q_flt_ps = - sin_theta_w_2.*v_f_2_al_flt_ps + cos_theta_w_2.*v_f_2_be_flt_ps;
    
    delta_theta_w_2 = atan2(v_f_2_q_flt_ps, v_f_2_d_flt_ps);
    delta_omega_w_2 = e_PLL_w_2 + PLL_Kp*delta_theta_w_2;

    % --- PQ control
    
    p_f_2 = v_f_2_a.*i_f_2_a + v_f_2_b.*i_f_2_b + v_f_2_c.*i_f_2_c;
    q_f_2 = 1/sqrt(3) * (v_f_2_a.*(i_f_2_c-i_f_2_b) + v_f_2_b.*(i_f_2_a-i_f_2_c) + v_f_2_c.*(i_f_2_b-i_f_2_a));

    i_f_2_d_ref = e_TPQ_2_d + TPQ_Kp * (p_f_2_ref - p_f_2);
    i_f_2_q_ref = e_TPQ_2_q - TPQ_Kp * (q_f_2_ref - q_f_2);

    i_f_2_al_ref = cos_theta_w_2.*i_f_2_d_ref - sin_theta_w_2.*i_f_2_q_ref;
    i_f_2_be_ref = sin_theta_w_2.*i_f_2_d_ref + cos_theta_w_2.*i_f_2_q_ref;

    % --- AC control

    v_w_2_al_ref = v_f_2_al_flt + e_TAC2_2_al + TAC_Kp_2 * (i_f_2_al_ref - i_f_2_al);
    v_w_2_be_ref = v_f_2_be_flt + e_TAC2_2_be + TAC_Kp_2 * (i_f_2_be_ref - i_f_2_be);

    % --- Modulation indices
    
    m_w_2_al_ref = v_w_2_al_ref ./ (v_d_w_2/2);
    m_w_2_be_ref = v_w_2_be_ref ./ (v_d_w_2/2);

    if is_delayed
        m_w_2_al = m_w_2_al_ref_dlyd;
        m_w_2_be = m_w_2_be_ref_dlyd;
    else
        m_w_2_al = m_w_2_al_ref;
        m_w_2_be = m_w_2_be_ref;
    end
    
    % --- voltage modulation

    v_w_2_al = v_d_w_2/2 .* m_w_2_al;
    v_w_2_be = v_d_w_2/2 .* m_w_2_be;
    
    % --- direct current
    i_d_w_2 = 3/4 * (m_w_2_al.*i_f_2_al + m_w_2_be.*i_f_2_be);

    % _____________________________________________________________________
    %                                                                 TLC_3
    % --- Clarke transformations
    
    i_f_3_a = i_f_3_al;
    i_f_3_b = 1/2 * (-i_f_3_al + sqrt(3)*i_f_3_be);
    i_f_3_c = 1/2 * (-i_f_3_al - sqrt(3)*i_f_3_be);
    
    i_x_3_a = i_x_3_al;
    i_x_3_b = 1/2 * (-i_x_3_al + sqrt(3)*i_x_3_be);
    i_x_3_c = 1/2 * (-i_x_3_al - sqrt(3)*i_x_3_be);
    
    v_f_3_a = v_f_3_al;
    v_f_3_b = 1/2 * (-v_f_3_al + sqrt(3)*v_f_3_be);
    v_f_3_c = 1/2 * (-v_f_3_al - sqrt(3)*v_f_3_be);
    
    % --- Transformer

    v_t_3_al = v_x_3_al / k_TLC;
    v_t_3_be = v_x_3_be / k_TLC;
    
    i_y_3_al = i_t_3_al / k_TLC;
    i_y_3_be = i_t_3_be / k_TLC;

    % --- Positive and negative sequence extraction
    
    v_f_3_al_flt_ps = 1/2 * ( v_f_3_al_flt      - v_f_3_be_flt_quad);
    v_f_3_be_flt_ps = 1/2 * ( v_f_3_al_flt_quad + v_f_3_be_flt     );
    v_f_3_al_flt_ns = 1/2 * ( v_f_3_al_flt      + v_f_3_be_flt_quad);
    v_f_3_be_flt_ns = 1/2 * (-v_f_3_al_flt_quad + v_f_3_be_flt     );
    
    % --- Park transformation and PLL
    
    cos_theta_eps_w_3 = cos(theta_eps_w_3);
    sin_theta_eps_w_3 = sin(theta_eps_w_3);

    cos_theta_w_3 = cos_omega_1_t.*cos_theta_eps_w_3 - sin_omega_1_t.*sin_theta_eps_w_3;
    sin_theta_w_3 = sin_omega_1_t.*cos_theta_eps_w_3 + cos_omega_1_t.*sin_theta_eps_w_3;

    % the PLL synchronises on the filtered positive-sequence voltage:
    v_f_3_d_flt_ps =   cos_theta_w_3.*v_f_3_al_flt_ps + sin_theta_w_3.*v_f_3_be_flt_ps;
    v_f_3_q_flt_ps = - sin_theta_w_3.*v_f_3_al_flt_ps + cos_theta_w_3.*v_f_3_be_flt_ps;
    
    delta_theta_w_3 = atan2(v_f_3_q_flt_ps, v_f_3_d_flt_ps);
    delta_omega_w_3 = e_PLL_w_3 + PLL_Kp*delta_theta_w_3;

    % --- PQ control
    
    p_f_3 = v_f_3_a.*i_f_3_a + v_f_3_b.*i_f_3_b + v_f_3_c.*i_f_3_c;
    q_f_3 = 1/sqrt(3) * (v_f_3_a.*(i_f_3_c-i_f_3_b) + v_f_3_b.*(i_f_3_a-i_f_3_c) + v_f_3_c.*(i_f_3_b-i_f_3_a));

    i_f_3_d_ref = e_TPQ_3_d + TPQ_Kp * (p_f_3_ref - p_f_3);
    i_f_3_q_ref = e_TPQ_3_q - TPQ_Kp * (q_f_3_ref - q_f_3);

    i_f_3_al_ref = cos_theta_w_3.*i_f_3_d_ref - sin_theta_w_3.*i_f_3_q_ref;
    i_f_3_be_ref = sin_theta_w_3.*i_f_3_d_ref + cos_theta_w_3.*i_f_3_q_ref;

    % --- AC control

    v_w_3_al_ref = v_f_3_al_flt + e_TAC2_3_al + TAC_Kp_3 * (i_f_3_al_ref - i_f_3_al);
    v_w_3_be_ref = v_f_3_be_flt + e_TAC2_3_be + TAC_Kp_3 * (i_f_3_be_ref - i_f_3_be);

    % --- Modulation indices
    
    m_w_3_al_ref = v_w_3_al_ref ./ (v_d_w_3/2);
    m_w_3_be_ref = v_w_3_be_ref ./ (v_d_w_3/2);

    if is_delayed
        m_w_3_al = m_w_3_al_ref_dlyd;
        m_w_3_be = m_w_3_be_ref_dlyd;
    else
        m_w_3_al = m_w_3_al_ref;
        m_w_3_be = m_w_3_be_ref;
    end

    % --- voltage modulation

    v_w_3_al = v_d_w_3/2 .* m_w_3_al;
    v_w_3_be = v_d_w_3/2 .* m_w_3_be;
    
    % --- direct current
    i_d_w_3 = 3/4 * (m_w_3_al.*i_f_3_al + m_w_3_be.*i_f_3_be);

    % _____________________________________________________________________
    %                                                          offshore PCC
    
    v_g_2_al = 2/3 * v_g_2_a - 1/3 * (v_g_2_b + v_g_2_c);
    v_g_2_be = 1/sqrt(3) * (v_g_2_b - v_g_2_c);

    i_x_sum_a = i_x_1_a + i_x_2_a + i_x_3_a;
    i_x_sum_b = i_x_1_b + i_x_2_b + i_x_3_b;
    i_x_sum_c = i_x_1_c + i_x_2_c + i_x_3_c;
    
end

if run_dxdt % defining dxdt = f(x,u)
    dxdt = [
        % _____________________________________________________________________
        %                                                                 TLC_1
        1/L_f_1 * (-R_f_1*i_f_1_al + v_w_1_al - v_f_1_al);              % i_f_1_al
        1/L_f_1 * (-R_f_1*i_f_1_be + v_w_1_be - v_f_1_be);              % i_f_1_be
        1/C_f_1 * (i_f_1_al - i_t_1_al);                                % v_f_1_al
        1/C_f_1 * (i_f_1_be - i_t_1_be);                                % v_f_1_be
        1/L_t_1 * (-R_t_1*i_t_1_al + v_f_1_al - v_t_1_al);              % i_t_1_al
        1/L_t_1 * (-R_t_1*i_t_1_be + v_f_1_be - v_t_1_be);              % i_t_1_be
        1/(C_x_1/2) * (i_y_1_al - i_x_1_al);                            % v_x_1_al
        1/(C_x_1/2) * (i_y_1_be - i_x_1_be);                            % v_x_1_be
        1/L_x_1 * (-R_x_1*i_x_1_al + v_x_1_al - v_g_2_al);              % i_x_1_al
        1/L_x_1 * (-R_x_1*i_x_1_be + v_x_1_be - v_g_2_be);              % i_x_1_be
        - omega_1 * v_f_1_al_flt;                                       % e_f_1_al_flt
        - omega_1 * v_f_1_be_flt;                                       % e_f_1_be_flt
        + omega_1 * (e_f_1_al_flt + v_f_1_al - v_f_1_al_flt);           % v_f_1_al_flt
        + omega_1 * (e_f_1_be_flt + v_f_1_be - v_f_1_be_flt);           % v_f_1_be_flt
        + omega_1 * (v_f_1_al - v_f_1_al_flt_quad);                     % e_f_1_al_flt_quad
        + omega_1 * (v_f_1_be - v_f_1_be_flt_quad);                     % e_f_1_be_flt_quad
        + omega_1 * (e_f_1_al_flt_quad - v_f_1_al_flt_quad);            % v_f_1_al_flt_quad
        + omega_1 * (e_f_1_be_flt_quad - v_f_1_be_flt_quad);            % v_f_1_be_flt_quad
        PLL_Ki * delta_theta_w_1;                                       % e_PLL_w_1
        delta_omega_w_1;                                                % theta_eps_w_1
        + TPQ_Ki * (p_f_1_ref - p_f_1);                                 % e_TPQ_1_d
        - TPQ_Ki * (q_f_1_ref - q_f_1);                                 % e_TPQ_1_q
        - omega_1 * e_TAC2_1_al;                                        % e_TAC_1_al
        - omega_1 * e_TAC2_1_be;                                        % e_TAC_1_be
        + omega_1 * e_TAC1_1_al + TAC_Kr_1 * (i_f_1_al_ref - i_f_1_al); % e_TAC_2_al
        + omega_1 * e_TAC1_1_be + TAC_Kr_1 * (i_f_1_be_ref - i_f_1_be); % e_TAC_2_be
        % _____________________________________________________________________
        %                                                                 TLC_2
        1/L_f_2 * (-R_f_2*i_f_2_al + v_w_2_al - v_f_2_al);              % i_f_2_al
        1/L_f_2 * (-R_f_2*i_f_2_be + v_w_2_be - v_f_2_be);              % i_f_2_be
        1/C_f_2 * (i_f_2_al - i_t_2_al);                                % v_f_2_al
        1/C_f_2 * (i_f_2_be - i_t_2_be);                                % v_f_2_be
        1/L_t_2 * (-R_t_2*i_t_2_al + v_f_2_al - v_t_2_al);              % i_t_2_al
        1/L_t_2 * (-R_t_2*i_t_2_be + v_f_2_be - v_t_2_be);              % i_t_2_be
        1/(C_x_2/2) * (i_y_2_al - i_x_2_al);                            % v_x_2_al
        1/(C_x_2/2) * (i_y_2_be - i_x_2_be);                            % v_x_2_be
        1/L_x_2 * (-R_x_2*i_x_2_al + v_x_2_al - v_g_2_al);              % i_x_2_al
        1/L_x_2 * (-R_x_2*i_x_2_be + v_x_2_be - v_g_2_be);              % i_x_2_be
        - omega_1 * v_f_2_al_flt;                                       % e_f_2_al_flt
        - omega_1 * v_f_2_be_flt;                                       % e_f_2_be_flt
        + omega_1 * (e_f_2_al_flt + v_f_2_al - v_f_2_al_flt);           % v_f_2_al_flt
        + omega_1 * (e_f_2_be_flt + v_f_2_be - v_f_2_be_flt);           % v_f_2_be_flt
        + omega_1 * (v_f_2_al - v_f_2_al_flt_quad);                     % e_f_2_al_flt_quad
        + omega_1 * (v_f_2_be - v_f_2_be_flt_quad);                     % e_f_2_be_flt_quad
        + omega_1 * (e_f_2_al_flt_quad - v_f_2_al_flt_quad);            % v_f_2_al_flt_quad
        + omega_1 * (e_f_2_be_flt_quad - v_f_2_be_flt_quad);            % v_f_2_be_flt_quad
        PLL_Ki * delta_theta_w_2;                                       % e_PLL_w_2
        delta_omega_w_2;                                                % theta_eps_w_2
        + TPQ_Ki * (p_f_2_ref - p_f_2);                                 % e_TPQ_2_d
        - TPQ_Ki * (q_f_2_ref - q_f_2);                                 % e_TPQ_2_q
        - omega_1 * e_TAC2_2_al;                                        % e_TAC_2_al
        - omega_1 * e_TAC2_2_be;                                        % e_TAC_2_be
        + omega_1 * e_TAC1_2_al + TAC_Kr_2 * (i_f_2_al_ref - i_f_2_al); % e_TAC_2_al
        + omega_1 * e_TAC1_2_be + TAC_Kr_2 * (i_f_2_be_ref - i_f_2_be); % e_TAC_2_be
        % _____________________________________________________________________
        %                                                                 TLC_3
        1/L_f_3 * (-R_f_3*i_f_3_al + v_w_3_al - v_f_3_al);              % i_f_3_al
        1/L_f_3 * (-R_f_3*i_f_3_be + v_w_3_be - v_f_3_be);              % i_f_3_be
        1/C_f_3 * (i_f_3_al - i_t_3_al);                                % v_f_3_al
        1/C_f_3 * (i_f_3_be - i_t_3_be);                                % v_f_3_be
        1/L_t_3 * (-R_t_3*i_t_3_al + v_f_3_al - v_t_3_al);              % i_t_3_al
        1/L_t_3 * (-R_t_3*i_t_3_be + v_f_3_be - v_t_3_be);              % i_t_3_be
        1/(C_x_3/2) * (i_y_3_al - i_x_3_al);                            % v_x_3_al
        1/(C_x_3/2) * (i_y_3_be - i_x_3_be);                            % v_x_3_be
        1/L_x_3 * (-R_x_3*i_x_3_al + v_x_3_al - v_g_2_al);              % i_x_3_al
        1/L_x_3 * (-R_x_3*i_x_3_be + v_x_3_be - v_g_2_be);              % i_x_3_be
        - omega_1 * v_f_3_al_flt;                                       % e_f_3_al_flt
        - omega_1 * v_f_3_be_flt;                                       % e_f_3_be_flt
        + omega_1 * (e_f_3_al_flt + v_f_3_al - v_f_3_al_flt);           % v_f_3_al_flt
        + omega_1 * (e_f_3_be_flt + v_f_3_be - v_f_3_be_flt);           % v_f_3_be_flt
        + omega_1 * (v_f_3_al - v_f_3_al_flt_quad);                     % e_f_3_al_flt_quad
        + omega_1 * (v_f_3_be - v_f_3_be_flt_quad);                     % e_f_3_be_flt_quad
        + omega_1 * (e_f_3_al_flt_quad - v_f_3_al_flt_quad);            % v_f_3_al_flt_quad
        + omega_1 * (e_f_3_be_flt_quad - v_f_3_be_flt_quad);            % v_f_3_be_flt_quad
        PLL_Ki * delta_theta_w_3;                                       % e_PLL_w_3
        delta_omega_w_3;                                                % theta_eps_w_3
        + TPQ_Ki * (p_f_3_ref - p_f_3);                                 % e_TPQ_3_d
        - TPQ_Ki * (q_f_3_ref - q_f_3);                                 % e_TPQ_3_q
        - omega_1 * e_TAC2_3_al;                                        % e_TAC_3_al
        - omega_1 * e_TAC2_3_be;                                        % e_TAC_3_be
        + omega_1 * e_TAC1_3_al + TAC_Kr_3 * (i_f_3_al_ref - i_f_3_al); % e_TAC_2_al
        + omega_1 * e_TAC1_3_be + TAC_Kr_3 * (i_f_3_be_ref - i_f_3_be); % e_TAC_2_be
        ];
end