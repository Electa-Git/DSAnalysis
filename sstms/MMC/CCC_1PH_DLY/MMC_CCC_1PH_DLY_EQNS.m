if get_values
    omega_1       = param.omega_1;
    R_a           = param.R_a;
    L_a           = param.L_a;
    C_a           = param.C_a;
    R_g           = param.R_g;
    R_e           = param.R_e;
    L_e           = param.L_e;
    AC_Kp         = param.AC_Kp;
    AC_Kr         = param.AC_Kr;
    CC_Kp         = param.CC_Kp;
    CC_Kr         = param.CC_Kr;

    v_g     = inputs.v_g;
    v_d     = inputs.v_d;
    i_s_ref = inputs.i_s_ref;
    
    i_s     = states.i_s;
    i_c     = states.i_c;
    v_c_u   = states.v_c_u;
    v_c_l   = states.v_c_l;
    i_d_flt = states.i_d_flt;
    eta_AC1 = states.eta_AC1;
    eta_AC2 = states.eta_AC2;
    eta_CC1 = states.eta_CC1;
    eta_CC2 = states.eta_CC2;

    n_u = delayed.n_u;
    n_l = delayed.n_l;
end

if run_algebra
    
    % ---  upper and lower-arm currents
    
    i_u = i_c + i_s/2;
    i_l = i_c - i_s/2;

    % ---  direct current
    
    i_d = i_c;
    
    % ---  instantaneous (single-phase) power and power losses
    
    p_g = v_g.*i_s;
    
    % --- alternating current control
    
    v_s_ref = v_g + eta_AC2 + AC_Kp * (i_s_ref - i_s);
    
    % ---  circulating current control
    
    i_c_ref = i_d_flt;
    v_c_ref = v_d/2 + eta_CC2 - CC_Kp * (i_c_ref - i_c);
    
    % --- modulation indices: uncompensated modulation

    n_u_ref = (-v_s_ref + v_c_ref) ./ v_d;
    n_l_ref = ( v_s_ref + v_c_ref) ./ v_d;
    
    if ~is_delayed % this close MUST be included
        n_u = n_u_ref;
        n_l = n_l_ref;
    end

    % --- inserted and driving voltages
    
    v_u = n_u .* v_c_u;
    v_l = n_l .* v_c_l;
    
    v_s = 1/2 * (-v_u + v_l);
    v_c = 1/2 * ( v_u + v_l);
end

if run_dxdt % defining dxdt = f(x,u)
    dxdt = [
        1/L_e * (-R_e*i_s + v_s - v_g);   % i_s
        1/L_a * (-R_a*i_c + v_d/2 - v_c); % i_c
        1/C_a * n_u .* i_u;               % v_c_u
        1/C_a * n_l .* i_l;               % v_c_l
        0.2*omega_1 * (i_d - i_d_flt);                  % i_d_flt
        -omega_1 * eta_AC2;                             % eta_AC1
        +omega_1 * eta_AC1 + AC_Kr * (i_s_ref - i_s);   % eta_AC2
        -2*omega_1 * eta_CC2;                           % eta_CC1
        +2*omega_1 * eta_CC1 - CC_Kr * (i_c_ref - i_c); % eta_CC2
        ];
end
