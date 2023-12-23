function inputs = MMC_CCC_1PH_DLY_uOft(t, DSA)
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
    I_s_d_ref = +0.8; % in pu
    
    % filling in numerical values at time t:
    inputs.v_g = V_g_ph_pk * sin(omega_1*t);
    inputs.v_d = V_d_pk * one;
    inputs.i_s_ref = I_s_d_ref * sin(omega_1*t);
end