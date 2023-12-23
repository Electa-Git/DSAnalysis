function DSA = RLC_DATA(DSA)
    
    % ---------- GENERAL NETWORK DATA ---------- %
    
    f_1     = 1;        % [Hz], fundamental frequency
    omega_1 = 2*pi*f_1; % [rad/s]
    S_nom   = 1;        % [VA]
    V_nom   = 1;        % [V] peak
    
    % ---------- BASE VALUES CALCULATION ---------- %
    
    S_base = S_nom;           % VA
    E_base = S_base;          % VA
    V_base = V_nom;           % V
    I_base = S_base / V_base; % A
    Z_base = V_base / I_base; % ohm
    
    base.S_base = S_base;
    base.E_base = E_base;
    base.V_base = V_base;
    base.I_base = I_base;
    base.Z_base = Z_base;
    
    % ---------- CIRCUIT PARAMETERS ---------- %
    
    R = 0.8;  % [ohm]
    L = 1;    % [H]
    C = 1e-3; % [F]
    
    % parameters are stored in structure "param", in per unit:
    param.R = R / Z_base;
    param.L = L / Z_base;
    param.C = C * Z_base;
    param.omega_1 = omega_1;
    
    % ---------- CONTROL PARAMETERS ---------- %
    
    % none
    
    % ---------- EXTERNAL HD-SOURCES ---------- %
    
    % none
    
    % ---------- INPUT PARAMETERS ---------- %
    
    V_s = V_base;
    
    % input parameters are also stored in structure "param", in per unit:
    param.V_s = V_s / V_base;
    
    % ---------- DELAY PARAMETERS ---------- %
    
    % none
    
    % ---------- DATA STRUCTURES ---------- %
    
    data.f_1   = f_1;
    data.base  = base;
    data.param = param;
    DSA.data   = data;
end