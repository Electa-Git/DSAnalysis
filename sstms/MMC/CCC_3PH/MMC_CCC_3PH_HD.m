function HDi = MMC_CCC_3PH_HD(DSA, HDi)
    % This function initialises the harmonic spectrum of:
    %   - the input   variables
    %   - the state   variables (initial guess for the solver)
    %   - the delayed variables (initial guess for the solver)
    % 
    % The initial spectra are stored in HDi, where 'i' refers to
    % initialisation.
    % 
    % -- Examples:
    %
    % 1- Defining a three-phase fundamental-frequency positive-sequence set
    % of variables:
    %   HDi.variable(h2i([-1 1])) = DSA_sine(amplitude, 0);
    %   HDi.variable(h2i([-1 1])) = DSA_sine(amplitude, -120);
    %   HDi.variable(h2i([-1 1])) = DSA_sine(amplitude, -240);
    % 
    %   (note that h_m must be at least 1 for this to work)
    %
    % 2- Defining a three-phase second harmonic negative-sequence set
    % of variables:
    %   HDi.variable(h2i([-2 2])) = DSA_sine(amplitude, 0);
    %   HDi.variable(h2i([-2 2])) = DSA_sine(amplitude, +120);
    %   HDi.variable(h2i([-2 2])) = DSA_sine(amplitude, +240);
    %
    %   (note that h_m must be at least 2 for this to work)
    % 
    % 3- Defining a constant variable:
    %   HDi.variable(h2i(0)) = amplitude;
    % 
    % -- Useful functions:
    % 
    %   h2i() takes harmonic indices and returns corresponding matlab 
    %   indices.
    %
    %   DSA_sine and DSA_cosine take peak amplitude and phase angle of a 
    %   single-frequency sinusoidal waveform and return negative and 
    %   positive frequency Fourier coefficients. The phase angle is 
    %   given in degrees.
    
    % extraction:
    h2i    = DSA.fcts.h2i_h_m;  % harmonic to matlab index transformation
    param  = DSA.data.param;    % parameters structure
    
    % _____________________________________________________________________
    %                                                  edit below this line
    
    % defining fundamental frequency of the current set of values:
    HDi.f_b = 50; % [Hz]
    
    % extracting parameters:
    V_g_ph_pk = param.V_g_ph_pk;
    V_d_pk    = param.V_d_pk;
    P_g_ref   = param.P_g_ref;
    Q_g_ref   = param.Q_g_ref;
    
    % --- Input variables:

    HDi.v_g_a(h2i([-1 1]))           = DSA_sine(V_g_ph_pk, 0);
    HDi.v_g_b(h2i([-1 1]))           = DSA_sine(V_g_ph_pk, -120);
    HDi.v_g_c(h2i([-1 1]))           = DSA_sine(V_g_ph_pk, -240);
    HDi.v_d(h2i(0))                  = V_d_pk;
    HDi.p_g_ref(h2i(0))              = P_g_ref;
    HDi.q_g_ref(h2i(0))              = Q_g_ref;
    HDi.cos_omega_1_t(h2i([-1 1]))   = DSA_cosine(1, 0);
    HDi.sin_omega_1_t(h2i([-1 1]))   = DSA_sine(1, 0);
    
    % --- State variables:
    % (just initial guesses)
    
    HDi.i_s_al(h2i([-1 1]))          = DSA_sine(1, 0);
    HDi.i_s_be(h2i([-1 1]))          = DSA_sine(1, 90);
    HDi.i_c_a(h2i(0))                = 1;
    HDi.i_c_b(h2i(0))                = 1;
    HDi.i_c_c(h2i(0))                = 1;
    HDi.v_c_u_a(h2i(0))              = 1;
    HDi.v_c_u_b(h2i(0))              = 1;
    HDi.v_c_u_c(h2i(0))              = 1;
    HDi.v_c_l_a(h2i(0))              = 1;
    HDi.v_c_l_b(h2i(0))              = 1;
    HDi.v_c_l_c(h2i(0))              = 1;
    HDi.i_d_flt(h2i(0))              = 1;
    HDi.v_g_al_flt(h2i([-1 1]))      = DSA_sine(1, 0);
    HDi.v_g_be_flt(h2i([-1 1]))      = DSA_sine(1, 90);
    HDi.v_g_al_flt_quad(h2i([-1 1])) = DSA_sine(1, 0);
    HDi.v_g_be_flt_quad(h2i([-1 1])) = DSA_sine(1, -90);
    
    % --- Delayed variables:
    % (just initial guesses)
    
    % none
end