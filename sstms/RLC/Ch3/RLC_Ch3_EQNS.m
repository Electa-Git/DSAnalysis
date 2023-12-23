if get_values
    omega_1 = param.omega_1;
    R       = param.R;
    L       = param.L;
    C       = param.C;

    v_s = inputs.v_s;
    i_l = states.i_l;
    v_c = states.v_c;
end

if run_algebra
    % none
end

if run_dxdt % defining dxdt = f(x,u)
    dxdt = [
        1/L * (-R*i_l - v_c + v_s); % i_l
        1/C * i_l;                  % v_c
        ];
end