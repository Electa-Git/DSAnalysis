function inputs = MMC_CCC_3PH_uOft(t, DSA)
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
    omega_1   = param.omega_1;
    V_g_ph_pk = param.V_g_ph_pk;
    V_d_pk    = param.V_d_pk;
    P_g_ref   = param.P_g_ref;
    Q_g_ref   = param.Q_g_ref;
    
    % filling in numerical values at time(s) t:
    inputs.v_g_a         = V_g_ph_pk * sin(omega_1*t - 0/3*pi);
    inputs.v_g_b         = V_g_ph_pk * sin(omega_1*t - 2/3*pi);
    inputs.v_g_c         = V_g_ph_pk * sin(omega_1*t - 4/3*pi);
    inputs.v_d           = V_d_pk * one;
    inputs.p_g_ref       = P_g_ref * one;
    inputs.q_g_ref       = Q_g_ref * one;
    inputs.cos_omega_1_t = cos(omega_1*t);
    inputs.sin_omega_1_t = sin(omega_1*t);
    
    % step types:
    % 0: do nothing
    % 1: adding negative-sequence AC voltage at fundamental frequency
    % 2: active power reference step
    
    step_type = 0;
    step_time = 0.02; % [s]

    mask = t >= step_time;
    if any(mask)
        switch step_type
            case 0
                % do nothing
            case 1
                inputs.v_g_a(mask) = V_g_ph_pk * sin(omega_1*t(mask) - 0/3*pi) + 0.02*V_g_ph_pk * sin(omega_1*t(mask) + 0/3*pi);
                inputs.v_g_b(mask) = V_g_ph_pk * sin(omega_1*t(mask) - 2/3*pi) + 0.02*V_g_ph_pk * sin(omega_1*t(mask) + 2/3*pi);
                inputs.v_g_c(mask) = V_g_ph_pk * sin(omega_1*t(mask) - 4/3*pi) + 0.02*V_g_ph_pk * sin(omega_1*t(mask) + 4/3*pi);
            case 2
                inputs.p_g_ref(mask) = 1.1*P_g_ref;
        end
    end
end