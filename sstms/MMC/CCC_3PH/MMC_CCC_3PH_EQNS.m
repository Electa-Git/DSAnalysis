if get_values
    omega_1 = param.omega_1;
    R_a     = param.R_a;
    L_a     = param.L_a;
    C_a     = param.C_a;
    R_e     = param.R_e;
    L_e     = param.L_e;
    AC_Kp   = param.AC_Kp;
    AC_Kr   = param.AC_Kr;
    CC_Kp   = param.CC_Kp;
    CC_Kr   = param.CC_Kr;
    PQ_Kp   = param.PQ_Kp;
    PQ_Ki   = param.PQ_Ki;
    PLL_Kp  = param.PLL_Kp;
    PLL_Ki  = param.PLL_Ki;
    
    v_g_a         = inputs.v_g_a;
    v_g_b         = inputs.v_g_b;
    v_g_c         = inputs.v_g_c;
    v_d           = inputs.v_d;
    p_g_ref       = inputs.p_g_ref;
    q_g_ref       = inputs.q_g_ref;
    cos_omega_1_t = inputs.cos_omega_1_t;
    sin_omega_1_t = inputs.sin_omega_1_t;
    
    i_s_al          = states.i_s_al;
    i_s_be          = states.i_s_be;
    i_c_a           = states.i_c_a;
    i_c_b           = states.i_c_b;
    i_c_c           = states.i_c_c;
    v_c_u_a         = states.v_c_u_a;
    v_c_u_b         = states.v_c_u_b;
    v_c_u_c         = states.v_c_u_c;
    v_c_l_a         = states.v_c_l_a;
    v_c_l_b         = states.v_c_l_b;
    v_c_l_c         = states.v_c_l_c;
    i_d_flt         = states.i_d_flt;
    eta_flt_al      = states.eta_flt_al;
    eta_flt_be      = states.eta_flt_be;
    v_g_al_flt      = states.v_g_al_flt;
    v_g_be_flt      = states.v_g_be_flt;
    eta_quad_al     = states.eta_quad_al;
    eta_quad_be     = states.eta_quad_be;
    v_g_al_flt_quad = states.v_g_al_flt_quad;
    v_g_be_flt_quad = states.v_g_be_flt_quad;
    eta_PLL         = states.eta_PLL;
    theta_eps       = states.theta_eps;
    eta_PQ_d        = states.eta_PQ_d;
    eta_PQ_q        = states.eta_PQ_q;
    eta_AC1_al      = states.eta_AC1_al;
    eta_AC1_be      = states.eta_AC1_be;
    eta_AC2_al      = states.eta_AC2_al;
    eta_AC2_be      = states.eta_AC2_be;
    eta_CC1_a       = states.eta_CC1_a;
    eta_CC1_b       = states.eta_CC1_b;
    eta_CC1_c       = states.eta_CC1_c;
    eta_CC2_a       = states.eta_CC2_a;
    eta_CC2_b       = states.eta_CC2_b;
    eta_CC2_c       = states.eta_CC2_c;
end

if run_algebra
    
    % --- Clarke transformations

    i_s_a = i_s_al;
    i_s_b = 1/2 * (-i_s_al + sqrt(3)*i_s_be);
    i_s_c = 1/2 * (-i_s_al - sqrt(3)*i_s_be);

    v_g_al = 2/3 * v_g_a - 1/3 * (v_g_b + v_g_c);
    v_g_be = 1/sqrt(3) * (v_g_b - v_g_c);
    
    % --- Current calculations

    i_u_a = i_c_a + i_s_a/2;
    i_u_b = i_c_b + i_s_b/2;
    i_u_c = i_c_c + i_s_c/2;
    i_l_a = i_c_a - i_s_a/2;
    i_l_b = i_c_b - i_s_b/2;
    i_l_c = i_c_c - i_s_c/2;

    i_d = i_c_a + i_c_b + i_c_c;
    
    % --- Positive and negative sequence extraction
    
    v_g_al_flt_ps = 1/2 * ( v_g_al_flt      - v_g_be_flt_quad);
    v_g_be_flt_ps = 1/2 * ( v_g_al_flt_quad + v_g_be_flt     );
    v_g_al_flt_ns = 1/2 * ( v_g_al_flt      + v_g_be_flt_quad);
    v_g_be_flt_ns = 1/2 * (-v_g_al_flt_quad + v_g_be_flt     );
    
    % --- Park transformation and PLL
    
    cos_theta_eps = cos(theta_eps);
    sin_theta_eps = sin(theta_eps);

    cos_theta = cos_omega_1_t.*cos_theta_eps - sin_omega_1_t.*sin_theta_eps;
    sin_theta = sin_omega_1_t.*cos_theta_eps + sin_theta_eps.*cos_omega_1_t;

    % the PLL synchronises on the filtered positive-sequence voltage:
    v_g_d_flt =   cos_theta.*v_g_al_flt_ps + sin_theta.*v_g_be_flt_ps;
    v_g_q_flt = - sin_theta.*v_g_al_flt_ps + cos_theta.*v_g_be_flt_ps;
    
    delta_theta = atan2(v_g_q_flt, v_g_d_flt);
    delta_omega = eta_PLL + PLL_Kp*delta_theta;
    
    % ---  PQ control

    p_g = v_g_a.*i_s_a + v_g_b.*i_s_b + v_g_c.*i_s_c;
    q_g = 1/sqrt(3) * (v_g_a.*(i_s_c-i_s_b) + v_g_b.*(i_s_a-i_s_c) + v_g_c.*(i_s_b-i_s_a));

    i_s_d_ref = eta_PQ_d + PQ_Kp * (p_g_ref - p_g);
    i_s_q_ref = eta_PQ_q - PQ_Kp * (q_g_ref - q_g);

    i_s_al_ref = cos_theta.*i_s_d_ref - sin_theta.*i_s_q_ref;
    i_s_be_ref = sin_theta.*i_s_d_ref + cos_theta.*i_s_q_ref;
    
    % --- AC control

    v_s_al_ref = v_g_al_flt + eta_AC2_al + AC_Kp * (i_s_al_ref - i_s_al);
    v_s_be_ref = v_g_be_flt + eta_AC2_be + AC_Kp * (i_s_be_ref - i_s_be);

    v_s_a_ref = v_s_al_ref;
    v_s_b_ref = 1/2 * (-v_s_al_ref + sqrt(3)*v_s_be_ref);
    v_s_c_ref = 1/2 * (-v_s_al_ref - sqrt(3)*v_s_be_ref);
    
    % --- CC control

    i_c_ref = i_d_flt;
    
    v_c_a_ref = v_d/2 + eta_CC2_a - CC_Kp * (i_c_ref - i_c_a);
    v_c_b_ref = v_d/2 + eta_CC2_b - CC_Kp * (i_c_ref - i_c_b);
    v_c_c_ref = v_d/2 + eta_CC2_c - CC_Kp * (i_c_ref - i_c_c);

    % --- modulation indices
    
    n_u_a_ref = (-v_s_a_ref + v_c_a_ref) ./ v_d;
    n_u_b_ref = (-v_s_b_ref + v_c_b_ref) ./ v_d;
    n_u_c_ref = (-v_s_c_ref + v_c_c_ref) ./ v_d;
    n_l_a_ref = ( v_s_a_ref + v_c_a_ref) ./ v_d;
    n_l_b_ref = ( v_s_b_ref + v_c_b_ref) ./ v_d;
    n_l_c_ref = ( v_s_c_ref + v_c_c_ref) ./ v_d;

    % no delays:
    n_u_a = n_u_a_ref;
    n_u_b = n_u_b_ref;
    n_u_c = n_u_c_ref;
    n_l_a = n_l_a_ref;
    n_l_b = n_l_b_ref;
    n_l_c = n_l_c_ref;
    
    % --- inserted and driving voltages
    
    v_u_a = n_u_a.*v_c_u_a;
    v_u_b = n_u_b.*v_c_u_b;
    v_u_c = n_u_c.*v_c_u_c;
    v_l_a = n_l_a.*v_c_l_a;
    v_l_b = n_l_b.*v_c_l_b;
    v_l_c = n_l_c.*v_c_l_c;
    
    v_s_a = 1/2 * (-v_u_a + v_l_a);
    v_s_b = 1/2 * (-v_u_b + v_l_b);
    v_s_c = 1/2 * (-v_u_c + v_l_c);
    v_c_a = 1/2 * ( v_u_a + v_l_a);
    v_c_b = 1/2 * ( v_u_b + v_l_b);
    v_c_c = 1/2 * ( v_u_c + v_l_c);
    
    v_s_al = 2/3 * v_s_a - 1/3 * (v_s_b + v_s_c);
    v_s_be = 1/sqrt(3) * (v_s_b - v_s_c);
    
    v_c_sig_a = (v_c_u_a + v_c_l_a)/2;
    v_c_sig_b = (v_c_u_b + v_c_l_b)/2;
    v_c_sig_c = (v_c_u_c + v_c_l_c)/2;
    v_c_del_a = (v_c_u_a - v_c_l_a)/2;
    v_c_del_b = (v_c_u_b - v_c_l_b)/2;
    v_c_del_c = (v_c_u_c - v_c_l_c)/2;
end

if run_dxdt % defining dxdt = f(x,u)
    dxdt = [
        1/L_e * (v_s_al - v_g_al - R_e*i_s_al);                 % i_s_al
        1/L_e * (v_s_be - v_g_be - R_e*i_s_be);                 % i_s_be
        1/L_a * (v_d/2 - v_c_a - R_a*i_c_a);                    % i_c_a
        1/L_a * (v_d/2 - v_c_b - R_a*i_c_b);                    % i_c_b
        1/L_a * (v_d/2 - v_c_c - R_a*i_c_c);                    % i_c_c
        1/C_a * n_u_a.*i_u_a;                                   % v_c_u_a
        1/C_a * n_u_b.*i_u_b;                                   % v_c_u_b
        1/C_a * n_u_c.*i_u_c;                                   % v_c_u_c
        1/C_a * n_l_a.*i_l_a;                                   % v_c_l_a
        1/C_a * n_l_b.*i_l_b;                                   % v_c_l_b
        1/C_a * n_l_c.*i_l_c;                                   % v_c_l_c
        
        +0.1*omega_1 * (i_d - i_d_flt);                         % i_d_flt
        -omega_1 * v_g_al_flt;                                  % eta_FLT_al
        -omega_1 * v_g_be_flt;                                  % eta_FLT_be
        +omega_1 * (eta_flt_al + v_g_al - v_g_al_flt);          % v_g_al_flt
        +omega_1 * (eta_flt_be + v_g_be - v_g_be_flt);          % v_g_be_flt
        +omega_1 * (v_g_al - v_g_al_flt_quad);                  % eta_QUAD_al
        +omega_1 * (v_g_be - v_g_be_flt_quad);                  % eta_QUAD_be
        +omega_1 * (eta_quad_al - v_g_al_flt_quad);             % v_g_al_flt_quad
        +omega_1 * (eta_quad_be - v_g_be_flt_quad);             % v_g_be_flt_quad
        
        PLL_Ki * delta_theta;                                   % eta_PLL
        delta_omega;                                            % theta_eps
        + PQ_Ki * (p_g_ref - p_g);                              % eta_PQ_d
        - PQ_Ki * (q_g_ref - q_g);                              % eta_PQ_q
        -omega_1 * eta_AC2_al;                                  % eta_AC1_al
        -omega_1 * eta_AC2_be;                                  % eta_AC1_be
        omega_1 * eta_AC1_al + AC_Kr * (i_s_al_ref - i_s_al);   % eta_AC2_al
        omega_1 * eta_AC1_be + AC_Kr * (i_s_be_ref - i_s_be);   % eta_AC2_be
        -2*omega_1 * eta_CC2_a;                                 % eta_CC1_a
        -2*omega_1 * eta_CC2_b;                                 % eta_CC1_b
        -2*omega_1 * eta_CC2_c;                                 % eta_CC1_c
        2*omega_1 * eta_CC1_a - CC_Kr * (i_c_ref - i_c_a);      % eta_CC2_a
        2*omega_1 * eta_CC1_b - CC_Kr * (i_c_ref - i_c_b);      % eta_CC2_b
        2*omega_1 * eta_CC1_c - CC_Kr * (i_c_ref - i_c_c);      % eta_CC2_c
        ];
end