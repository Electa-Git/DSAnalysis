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
    k_TLC    = param.k_TLC;
    PLL_Kp   = param.PLL_Kp;
    PLL_Ki   = param.PLL_Ki;
    TAC_Kp_1 = param.TAC_Kp_1;
    TAC_Kr_1 = param.TAC_Kr_1;
    TPQ_Kp   = param.TPQ_Kp;
    TPQ_Ki   = param.TPQ_Ki;

    v_g_2_a       = inputs.v_g_2_a;
    v_g_2_b       = inputs.v_g_2_b;
    v_g_2_c       = inputs.v_g_2_c;
    p_f_1_ref     = inputs.p_f_1_ref;
    q_f_1_ref     = inputs.q_f_1_ref;
    v_d_w_1       = inputs.v_d_w_1;
    cos_omega_1_t = inputs.cos_omega_1_t;
    sin_omega_1_t = inputs.sin_omega_1_t;
    
    i_f_1_al          = states.i_f_1_al;
    i_f_1_be          = states.i_f_1_be;
    v_f_1_al          = states.v_f_1_al;
    v_f_1_be          = states.v_f_1_be;
    i_t_1_al          = states.i_t_1_al;
    i_t_1_be          = states.i_t_1_be;
    v_x_1_al          = states.v_x_1_al;
    v_x_1_be          = states.v_x_1_be;
    i_x_1_al          = states.i_x_1_al;
    i_x_1_be          = states.i_x_1_be;
    e_f_1_al_flt      = states.e_f_1_al_flt;
    e_f_1_be_flt      = states.e_f_1_be_flt;
    v_f_1_al_flt      = states.v_f_1_al_flt;
    v_f_1_be_flt      = states.v_f_1_be_flt;
    e_f_1_al_flt_quad = states.e_f_1_al_flt_quad;
    e_f_1_be_flt_quad = states.e_f_1_be_flt_quad;
    v_f_1_al_flt_quad = states.v_f_1_al_flt_quad;
    v_f_1_be_flt_quad = states.v_f_1_be_flt_quad;
    e_PLL_w_1         = states.e_PLL_w_1;
    theta_eps_w_1     = states.theta_eps_w_1;
    e_TPQ_1_d         = states.e_TPQ_1_d;
    e_TPQ_1_q         = states.e_TPQ_1_q;
    e_TAC1_1_al       = states.e_TAC1_1_al;
    e_TAC1_1_be       = states.e_TAC1_1_be;
    e_TAC2_1_al       = states.e_TAC2_1_al;
    e_TAC2_1_be       = states.e_TAC2_1_be;
    
    m_w_1_al_ref_dlyd = delayed.m_w_1_al_ref_dlyd;
    m_w_1_be_ref_dlyd = delayed.m_w_1_be_ref_dlyd;
end

if run_algebra

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

    v_g_2_al = 2/3 * v_g_2_a - 1/3 * (v_g_2_b + v_g_2_c);
    v_g_2_be = 1/sqrt(3) * (v_g_2_b - v_g_2_c);
    
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
end

if run_dxdt % defining dxdt = f(x,u)
    dxdt = [
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
        - omega_1 * e_TAC2_1_al;                                        % e_TAC1_1_al
        - omega_1 * e_TAC2_1_be;                                        % e_TAC1_1_be
        + omega_1 * e_TAC1_1_al + TAC_Kr_1 * (i_f_1_al_ref - i_f_1_al); % e_TAC2_1_al
        + omega_1 * e_TAC1_1_be + TAC_Kr_1 * (i_f_1_be_ref - i_f_1_be); % e_TAC2_1_be
        ];
end