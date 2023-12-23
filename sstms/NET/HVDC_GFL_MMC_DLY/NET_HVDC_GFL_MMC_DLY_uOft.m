function inputs = NET_HVDC_GFL_MMC_DLY_uOft(t, DSA)
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
    V_g_1_ph_pk = param.V_g_1_ph_pk;
    V_g_2_ph_pk = param.V_g_2_ph_pk;
    V_d_pk      = param.V_d_pk;
    Q_g_1_ref   = param.Q_g_1_ref;
    P_g_2_ref   = param.P_g_2_ref;
    Q_g_2_ref   = param.Q_g_2_ref;
    
    % filling in numerical values at time t:
    inputs.v_g_1_a       = V_g_1_ph_pk * sin(omega_1*t - 0/3*pi);
    inputs.v_g_1_b       = V_g_1_ph_pk * sin(omega_1*t - 2/3*pi);
    inputs.v_g_1_c       = V_g_1_ph_pk * sin(omega_1*t - 4/3*pi);
    inputs.v_g_2_a       = V_g_2_ph_pk * sin(omega_1*t - 0/3*pi);
    inputs.v_g_2_b       = V_g_2_ph_pk * sin(omega_1*t - 2/3*pi);
    inputs.v_g_2_c       = V_g_2_ph_pk * sin(omega_1*t - 4/3*pi);
    inputs.v_d_1_ref     = V_d_pk * one;
    inputs.q_g_1_ref     = Q_g_1_ref * one;
    inputs.p_g_2_ref     = P_g_2_ref * one;
    inputs.q_g_2_ref     = Q_g_2_ref * one;
    inputs.cos_omega_1_t = cos(omega_1*t);
    inputs.sin_omega_1_t = sin(omega_1*t);
    
    % step of external source or control reference:
    % --1 choose step type
    % --2 choose step time
    % --3 check the implementation of the step: masks are used to overwrite
    % all time values occuring after the step time

    % step types:
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
                inputs.v_g_1_a(mask) = V_g_1_ph_pk * sin(omega_1*t(mask) - 0/3*pi) + 0.02*V_g_1_ph_pk * sin(omega_1*t(mask) + 0/3*pi); % v_g_1_a
                inputs.v_g_1_b(mask) = V_g_1_ph_pk * sin(omega_1*t(mask) - 2/3*pi) + 0.02*V_g_1_ph_pk * sin(omega_1*t(mask) + 2/3*pi); % v_g_1_b
                inputs.v_g_1_c(mask) = V_g_1_ph_pk * sin(omega_1*t(mask) - 4/3*pi) + 0.02*V_g_1_ph_pk * sin(omega_1*t(mask) + 4/3*pi); % v_g_1_c
            case 2
                inputs.p_g_2_ref(mask) = 1.1*P_g_2_ref; % p_g_2_ref
            otherwise
                error('undefined step type')
        end
    end
end