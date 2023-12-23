function inputs = RLC_Ch3_uOft(t, DSA)
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
    omega_1 = param.omega_1;
    V_s     = param.V_s;
    
    % filling in numerical values at time t:
    inputs.v_s = V_s + 0.5*V_s*sin(omega_1*t); % v_s

    % adding a disturbance:
    step_time = +inf; % [s]
    
    mask = t >= step_time;
    if any(mask)
        inputs.v_s(mask) = V_s + 0.5*V_s*sin(omega_1*t(mask)) + 0.05*V_s*sin(5*omega_1*t(mask)); % v_s
    end
end