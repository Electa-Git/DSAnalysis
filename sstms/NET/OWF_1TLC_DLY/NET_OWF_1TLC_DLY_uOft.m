function inputs = NET_OWF_1TLC_DLY_uOft(t, DSA)
    % This function returns the value of input variables at time t.
    % 
    % Time t is either a scalar or a vector.
    % If t is a scalar, inputs will be scalars.
    % If t is a vector, inputs will be vectors.
    %
    % Steps and disturbances can be defined in this file. The use of a time
    % mask is then necessary, as it enables working with t as either a 
    % scalar or a vector.
    
    % parameters structure:
    param = DSA.data.param;
    
    % vector of ones for constant inputs:
    one = ones(size(t));
    
    % _____________________________________________________________________
    %                                                  edit below this line
    
    % getting parameter values:
    omega_1     = param.omega_1;
    V_g_2_ph_pk = param.V_g_2_ph_pk;
    V_d_w_pk    = param.V_d_w_pk;
    P_f_1_ref   = param.P_f_1_ref;
    Q_f_1_ref   = param.Q_f_1_ref;
    
    % filling in numerical values at time(s) t:
    inputs.v_g_2_a       = V_g_2_ph_pk * sin(omega_1*t - 0/3*pi);
    inputs.v_g_2_b       = V_g_2_ph_pk * sin(omega_1*t - 2/3*pi);
    inputs.v_g_2_c       = V_g_2_ph_pk * sin(omega_1*t - 4/3*pi);
    inputs.p_f_1_ref     = P_f_1_ref * one;
    inputs.q_f_1_ref     = Q_f_1_ref * one;
    inputs.v_d_w_1       = V_d_w_pk * one;
    inputs.cos_omega_1_t = cos(omega_1*t);
    inputs.sin_omega_1_t = sin(omega_1*t);
end








