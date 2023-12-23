function HDi = NET_GFM_MMC_OWF_3TLC_DLY_HD(DSA, HDi)
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
    
    % getting parameter values:
    I_a_ph_pk   = param.I_a_ph_pk;
    V_d_pk      = param.V_d_pk;
    I_d_pk      = param.I_d_pk;
    V_g_2_ph_pk = param.V_g_2_ph_pk;
    I_g_2_ph_pk = param.I_g_2_ph_pk;
    V_w_ph_pk   = param.V_w_ph_pk;
    I_w_1_ph_pk = param.I_w_1_ph_pk;
    I_w_2_ph_pk = param.I_w_2_ph_pk;
    I_w_3_ph_pk = param.I_w_3_ph_pk;
    V_d_w_pk    = param.V_d_w_pk;
    P_f_1_ref   = param.P_f_1_ref;
    Q_f_1_ref   = param.Q_f_1_ref;
    P_f_2_ref   = param.P_f_2_ref;
    Q_f_2_ref   = param.Q_f_2_ref;
    P_f_3_ref   = param.P_f_3_ref;
    Q_f_3_ref   = param.Q_f_3_ref;
    
    % --- Input variables:
    
    HDi.v_d_2(h2i(0))              = V_d_pk;
    HDi.v_g_2_a_ref(h2i([-1 1]))   = DSA_sine(V_g_2_ph_pk, 0);
    HDi.v_g_2_b_ref(h2i([-1 1]))   = DSA_sine(V_g_2_ph_pk, -120);
    HDi.v_g_2_c_ref(h2i([-1 1]))   = DSA_sine(V_g_2_ph_pk, -240);
    HDi.p_f_1_ref(h2i(0))          = P_f_1_ref;
    HDi.q_f_1_ref(h2i(0))          = Q_f_1_ref;
    HDi.p_f_2_ref(h2i(0))          = P_f_2_ref;
    HDi.q_f_2_ref(h2i(0))          = Q_f_2_ref;
    HDi.p_f_3_ref(h2i(0))          = P_f_3_ref;
    HDi.q_f_3_ref(h2i(0))          = Q_f_3_ref;
    HDi.v_d_w_1(h2i(0))            = V_d_w_pk;
    HDi.v_d_w_2(h2i(0))            = V_d_w_pk;
    HDi.v_d_w_3(h2i(0))            = V_d_w_pk;
    HDi.cos_omega_1_t(h2i([-1 1])) = DSA_cosine(1, 0);
    HDi.sin_omega_1_t(h2i([-1 1])) = DSA_sine(1, 0);

    % --- State variables:
    % (just initial guesses)
    
    HDi.i_s_2_al(h2i([-1 1]))          = DSA_sine(I_a_ph_pk, 0);
    HDi.i_s_2_be(h2i([-1 1]))          = DSA_sine(I_a_ph_pk, -90);
    HDi.i_c_2_a(h2i(0))                = 1;
    HDi.i_c_2_b(h2i(0))                = 1;
    HDi.i_c_2_c(h2i(0))                = 1;
    HDi.v_c_u_2_a(h2i(0))              = V_d_pk;
    HDi.v_c_u_2_b(h2i(0))              = V_d_pk;
    HDi.v_c_u_2_c(h2i(0))              = V_d_pk;
    HDi.v_c_l_2_a(h2i(0))              = V_d_pk;
    HDi.v_c_l_2_b(h2i(0))              = V_d_pk;
    HDi.v_c_l_2_c(h2i(0))              = V_d_pk;
    HDi.e_g_2_al_flt(h2i([-1 1]))      = DSA_sine(V_g_2_ph_pk, 90);
    HDi.e_g_2_be_flt(h2i([-1 1]))      = DSA_sine(V_g_2_ph_pk, 0);
    HDi.v_g_2_al_flt(h2i([-1 1]))      = DSA_sine(V_g_2_ph_pk, 0);
    HDi.v_g_2_be_flt(h2i([-1 1]))      = DSA_sine(V_g_2_ph_pk, -90);
    HDi.e_g_2_al_flt_quad(h2i([-1 1])) = DSA_sine(V_g_2_ph_pk, 0);
    HDi.e_g_2_be_flt_quad(h2i([-1 1])) = DSA_sine(V_g_2_ph_pk, -90);
    HDi.v_g_2_al_flt_quad(h2i([-1 1])) = DSA_sine(V_g_2_ph_pk, -90);
    HDi.v_g_2_be_flt_quad(h2i([-1 1])) = DSA_sine(V_g_2_ph_pk, -180);
    HDi.e_MAV1_al(h2i([-1 1]))         = DSA_sine(0.01, 0);
    HDi.e_MAV1_be(h2i([-1 1]))         = DSA_sine(0.01, -90);
    HDi.e_MAV2_al(h2i([-1 1]))         = DSA_sine(0.01, 0);
    HDi.e_MAV2_be(h2i([-1 1]))         = DSA_sine(0.01, -90);
    HDi.e_MAC1_2_al(h2i([-1 1]))       = DSA_sine(0.01, 0);
    HDi.e_MAC1_2_be(h2i([-1 1]))       = DSA_sine(0.01, -90);
    HDi.e_MAC2_2_al(h2i([-1 1]))       = DSA_sine(0.01, 0);
    HDi.e_MAC2_2_be(h2i([-1 1]))       = DSA_sine(0.01, -90);
    HDi.e_MCC1_2_a(h2i(0))             = 0.01;
    HDi.e_MCC1_2_b(h2i(0))             = 0.01;
    HDi.e_MCC1_2_c(h2i(0))             = 0.01;
    HDi.e_MCC2_2_a(h2i(0))             = 0.01;
    HDi.e_MCC2_2_b(h2i(0))             = 0.01;
    HDi.e_MCC2_2_c(h2i(0))             = 0.01;
    
    HDi.v_g_2_a(h2i([-1 1]))           = DSA_sine(V_g_2_ph_pk, 0);
    HDi.v_g_2_b(h2i([-1 1]))           = DSA_sine(V_g_2_ph_pk, -120);
    HDi.v_g_2_c(h2i([-1 1]))           = DSA_sine(V_g_2_ph_pk, -240);

    % -- TCL 1 and collection cable 1
    HDi.i_f_1_al(h2i([-1 1]))          = DSA_sine(I_w_1_ph_pk, 0);
    HDi.i_f_1_be(h2i([-1 1]))          = DSA_sine(I_w_1_ph_pk, -90);
    HDi.v_f_1_al(h2i([-1 1]))          = DSA_sine(V_w_ph_pk, 0);
    HDi.v_f_1_be(h2i([-1 1]))          = DSA_sine(V_w_ph_pk, -90);
    HDi.i_t_1_al(h2i([-1 1]))          = DSA_sine(I_w_1_ph_pk, 0);
    HDi.i_t_1_be(h2i([-1 1]))          = DSA_sine(I_w_1_ph_pk, -90);
    HDi.v_x_1_al(h2i([-1 1]))          = DSA_sine(V_g_2_ph_pk, 0);
    HDi.v_x_1_be(h2i([-1 1]))          = DSA_sine(V_g_2_ph_pk, -90);
    HDi.i_x_1_al(h2i([-1 1]))          = DSA_sine(I_g_2_ph_pk, 0);
    HDi.i_x_1_be(h2i([-1 1]))          = DSA_sine(I_g_2_ph_pk, -90);
    HDi.e_f_1_al_flt(h2i([-1 1]))      = DSA_sine(V_w_ph_pk, 90);
    HDi.e_f_1_be_flt(h2i([-1 1]))      = DSA_sine(V_w_ph_pk, 0);
    HDi.v_f_1_al_flt(h2i([-1 1]))      = DSA_sine(V_w_ph_pk, 0);
    HDi.v_f_1_be_flt(h2i([-1 1]))      = DSA_sine(V_w_ph_pk, -90);
    HDi.e_f_1_al_flt_quad(h2i([-1 1])) = DSA_sine(V_w_ph_pk, 0);
    HDi.e_f_1_be_flt_quad(h2i([-1 1])) = DSA_sine(V_w_ph_pk, -90);
    HDi.v_f_1_al_flt_quad(h2i([-1 1])) = DSA_sine(V_w_ph_pk, -90);
    HDi.v_f_1_be_flt_quad(h2i([-1 1])) = DSA_sine(V_w_ph_pk, -180);
    HDi.e_PLL_w_1(h2i(0))              = 0.01;
    HDi.theta_eps_w_1(h2i(0))          = 0.01;
    HDi.e_TPQ_1_d(h2i(0))              = 0.01;
    HDi.e_TPQ_1_q(h2i(0))              = 0.01;
    HDi.e_TAC1_1_al(h2i([-1 1]))       = DSA_sine(0.01, 0);
    HDi.e_TAC1_1_be(h2i([-1 1]))       = DSA_sine(0.01, -90);
    HDi.e_TAC2_1_al(h2i([-1 1]))       = DSA_sine(0.01, 0);
    HDi.e_TAC2_1_be(h2i([-1 1]))       = DSA_sine(0.01, -90);
    % -- TCL 2 and collection cable 2
    HDi.i_f_2_al(h2i([-1 1]))          = DSA_sine(I_w_2_ph_pk, 0);
    HDi.i_f_2_be(h2i([-1 1]))          = DSA_sine(I_w_2_ph_pk, -90);
    HDi.v_f_2_al(h2i([-1 1]))          = DSA_sine(V_w_ph_pk, 0);
    HDi.v_f_2_be(h2i([-1 1]))          = DSA_sine(V_w_ph_pk, -90);
    HDi.i_t_2_al(h2i([-1 1]))          = DSA_sine(I_w_2_ph_pk, 0);
    HDi.i_t_2_be(h2i([-1 1]))          = DSA_sine(I_w_2_ph_pk, -90);
    HDi.v_x_2_al(h2i([-1 1]))          = DSA_sine(V_g_2_ph_pk, 0);
    HDi.v_x_2_be(h2i([-1 1]))          = DSA_sine(V_g_2_ph_pk, -90);
    HDi.i_x_2_al(h2i([-1 1]))          = DSA_sine(I_g_2_ph_pk, 0);
    HDi.i_x_2_be(h2i([-1 1]))          = DSA_sine(I_g_2_ph_pk, -90);
    HDi.e_f_2_al_flt(h2i([-1 1]))      = DSA_sine(V_w_ph_pk, 90);
    HDi.e_f_2_be_flt(h2i([-1 1]))      = DSA_sine(V_w_ph_pk, 0);
    HDi.v_f_2_al_flt(h2i([-1 1]))      = DSA_sine(V_w_ph_pk, 0);
    HDi.v_f_2_be_flt(h2i([-1 1]))      = DSA_sine(V_w_ph_pk, -90);
    HDi.e_f_2_al_flt_quad(h2i([-1 1])) = DSA_sine(V_w_ph_pk, 0);
    HDi.e_f_2_be_flt_quad(h2i([-1 1])) = DSA_sine(V_w_ph_pk, -90);
    HDi.v_f_2_al_flt_quad(h2i([-1 1])) = DSA_sine(V_w_ph_pk, -90);
    HDi.v_f_2_be_flt_quad(h2i([-1 1])) = DSA_sine(V_w_ph_pk, -180);
    HDi.e_PLL_w_2(h2i(0))              = 0.01;
    HDi.theta_eps_w_2(h2i(0))          = 0.01;
    HDi.e_TPQ_2_d(h2i(0))              = 0.01;
    HDi.e_TPQ_2_q(h2i(0))              = 0.01;
    HDi.e_TAC1_2_al(h2i([-1 1]))       = DSA_sine(0.01, 0);
    HDi.e_TAC1_2_be(h2i([-1 1]))       = DSA_sine(0.01, -90);
    HDi.e_TAC2_2_al(h2i([-1 1]))       = DSA_sine(0.01, 0);
    HDi.e_TAC2_2_be(h2i([-1 1]))       = DSA_sine(0.01, -90);
    % -- TCL 3 and collection cable 3
    HDi.i_f_3_al(h2i([-1 1]))          = DSA_sine(I_w_3_ph_pk, 0);
    HDi.i_f_3_be(h2i([-1 1]))          = DSA_sine(I_w_3_ph_pk, -90);
    HDi.v_f_3_al(h2i([-1 1]))          = DSA_sine(V_w_ph_pk, 0);
    HDi.v_f_3_be(h2i([-1 1]))          = DSA_sine(V_w_ph_pk, -90);
    HDi.i_t_3_al(h2i([-1 1]))          = DSA_sine(I_w_3_ph_pk, 0);
    HDi.i_t_3_be(h2i([-1 1]))          = DSA_sine(I_w_3_ph_pk, -90);
    HDi.v_x_3_al(h2i([-1 1]))          = DSA_sine(V_g_2_ph_pk, 0);
    HDi.v_x_3_be(h2i([-1 1]))          = DSA_sine(V_g_2_ph_pk, -90);
    HDi.i_x_3_al(h2i([-1 1]))          = DSA_sine(I_g_2_ph_pk, 0);
    HDi.i_x_3_be(h2i([-1 1]))          = DSA_sine(I_g_2_ph_pk, -90);
    HDi.e_f_3_al_flt(h2i([-1 1]))      = DSA_sine(V_w_ph_pk, 90);
    HDi.e_f_3_be_flt(h2i([-1 1]))      = DSA_sine(V_w_ph_pk, 0);
    HDi.v_f_3_al_flt(h2i([-1 1]))      = DSA_sine(V_w_ph_pk, 0);
    HDi.v_f_3_be_flt(h2i([-1 1]))      = DSA_sine(V_w_ph_pk, -90);
    HDi.e_f_3_al_flt_quad(h2i([-1 1])) = DSA_sine(V_w_ph_pk, 0);
    HDi.e_f_3_be_flt_quad(h2i([-1 1])) = DSA_sine(V_w_ph_pk, -90);
    HDi.v_f_3_al_flt_quad(h2i([-1 1])) = DSA_sine(V_w_ph_pk, -90);
    HDi.v_f_3_be_flt_quad(h2i([-1 1])) = DSA_sine(V_w_ph_pk, -180);
    HDi.e_PLL_w_3(h2i(0))              = 0.01;
    HDi.theta_eps_w_3(h2i(0))          = 0.01;
    HDi.e_TPQ_3_d(h2i(0))              = 0.01;
    HDi.e_TPQ_3_q(h2i(0))              = 0.01;
    HDi.e_TAC1_3_al(h2i([-1 1]))       = DSA_sine(0.01, 0);
    HDi.e_TAC1_3_be(h2i([-1 1]))       = DSA_sine(0.01, -90);
    HDi.e_TAC2_3_al(h2i([-1 1]))       = DSA_sine(0.01, 0);
    HDi.e_TAC2_3_be(h2i([-1 1]))       = DSA_sine(0.01, -90);
    
    % --- Delayed variables:
    % (just initial guesses)
    
    HDi.n_u_2_a_ref_dlyd(h2i(0))      = 0.5;
    HDi.n_u_2_b_ref_dlyd(h2i(0))      = 0.5;
    HDi.n_u_2_c_ref_dlyd(h2i(0))      = 0.5;
    HDi.n_l_2_a_ref_dlyd(h2i(0))      = 0.5;
    HDi.n_l_2_b_ref_dlyd(h2i(0))      = 0.5;
    HDi.n_l_2_c_ref_dlyd(h2i(0))      = 0.5;
    HDi.n_u_2_a_ref_dlyd(h2i([-1 1])) = DSA_sine(0.5,     -180);
    HDi.n_u_2_b_ref_dlyd(h2i([-1 1])) = DSA_sine(0.5, -120-180);
    HDi.n_u_2_c_ref_dlyd(h2i([-1 1])) = DSA_sine(0.5, -240-180);
    HDi.n_l_2_a_ref_dlyd(h2i([-1 1])) = DSA_sine(0.5,    0);
    HDi.n_l_2_b_ref_dlyd(h2i([-1 1])) = DSA_sine(0.5, -120);
    HDi.n_l_2_c_ref_dlyd(h2i([-1 1])) = DSA_sine(0.5, -240);

    HDi.m_w_1_al_ref_dlyd(h2i([-1 1])) = DSA_sine(1, 0);
    HDi.m_w_1_be_ref_dlyd(h2i([-1 1])) = DSA_sine(1, -90);
    HDi.m_w_2_al_ref_dlyd(h2i([-1 1])) = DSA_sine(1, 0);
    HDi.m_w_2_be_ref_dlyd(h2i([-1 1])) = DSA_sine(1, -90);
    HDi.m_w_3_al_ref_dlyd(h2i([-1 1])) = DSA_sine(1, 0);
    HDi.m_w_3_be_ref_dlyd(h2i([-1 1])) = DSA_sine(1, -90);
end