function DSA = MMC_CCC_3PH_DLY_VARS(DSA)
    % This function defines lists of variables related to the system.
    % 
    % Variables are defined for:
    %     - the large signal (lasi)  nonlinear system 
    %     - the small signal (smsi) linearised system
    % 
    % For lasi systems:
    %     - Names of inputs, states, outputs are provided by the user.
    %
    % For systems with pure delays, additional lists must be populated:
    %     - Names of variables that will be delayed ("todelay")
    %     - Names of corresponding delayed variables ("delayed")
    %     - Otherwise, leave these lists empty.
    %
    % For smsi systems:
    %     - Names of inputs and outputs are provided by the user.
    %
    % In the lasi sets of variables, the outputs can be any intermediate
    % variable that the user would like to display or retrieve from the 
    % calculations.
    % 
    % In the smsi sets of variables, the inputs and outputs define the
    % variables between which frequency responses can be calculated.
    % 
    % The order in which the variables are listed is important:
    %     - The state variables must be in the same order as in the dxdt
    %      vector in the EQNS file.
    %     - The delayed and todelay lists must have matching order.
    % 
    % The only difference between lasi and smsi systems are the
    % lists of inputs and outputs. The states, todelay and delayed
    % variables remain the same.
    % _____________________________________________________________________
    %                                                  edit below this line
    
    lasi_inputs = {
                'v_g_a';
                'v_g_b';
                'v_g_c';
                'v_d';
                'p_g_ref';
                'q_g_ref';
                'cos_omega_1_t';
                'sin_omega_1_t';
                };
    lasi_states = {
                'i_s_al';
                'i_s_be';
                'i_c_a';
                'i_c_b';
                'i_c_c';
                'v_c_u_a';
                'v_c_u_b';
                'v_c_u_c';
                'v_c_l_a';
                'v_c_l_b';
                'v_c_l_c';
                'i_d_flt';
                'eta_flt_al';
                'eta_flt_be';
                'v_g_al_flt';
                'v_g_be_flt';
                'eta_quad_al';
                'eta_quad_be';
                'v_g_al_flt_quad';
                'v_g_be_flt_quad';
                'eta_PLL';
                'theta_eps';
                'eta_PQ_d';
                'eta_PQ_q';
                'eta_AC1_al';
                'eta_AC1_be';
                'eta_AC2_al';
                'eta_AC2_be';
                'eta_CC1_a';
                'eta_CC1_b';
                'eta_CC1_c';
                'eta_CC2_a';
                'eta_CC2_b';
                'eta_CC2_c';
                };
    lasi_outputs = {
                'i_s_a';
                'i_s_b';
                'i_s_c';
                'i_u_a';
                'i_u_b';
                'i_u_c';
                'i_l_a';
                'i_l_b';
                'i_l_c';
                'i_s_d_ref';
                'i_s_q_ref';
                'v_g_al';
                'v_g_be';
                'v_g_al_flt_ps';
                'v_g_be_flt_ps';
                'v_g_al_flt_ns';
                'v_g_be_flt_ns';
                'v_g_d_flt';
                'v_g_q_flt';
                'cos_theta';
                'sin_theta';
                'cos_theta_eps';
                'sin_theta_eps';
                'delta_theta';
                'delta_omega';
                'i_s_al_ref';
                'i_s_be_ref';
                'p_g';
                'q_g';
                'i_d';
                'i_c_ref';
                'v_s_al_ref';
                'v_s_be_ref';
                'v_s_a_ref';
                'v_s_b_ref';
                'v_s_c_ref';
                'v_c_a_ref';
                'v_c_b_ref';
                'v_c_c_ref';
                'n_u_a_ref';
                'n_u_b_ref';
                'n_u_c_ref';
                'n_l_a_ref';
                'n_l_b_ref';
                'n_l_c_ref';
                'n_u_a';
                'n_u_b';
                'n_u_c';
                'n_l_a';
                'n_l_b';
                'n_l_c';
                'v_u_a';
                'v_u_b';
                'v_u_c';
                'v_l_a';
                'v_l_b';
                'v_l_c';
                'v_s_a';
                'v_s_b';
                'v_s_c';
                'v_c_a';
                'v_c_b';
                'v_c_c';
                'v_s_al';
                'v_s_be';
                'v_c_sig_a';
                'v_c_sig_b';
                'v_c_sig_c';
                'v_c_del_a';
                'v_c_del_b';
                'v_c_del_c';
                };
    lasi_todelay = {
                'n_u_a_ref';
                'n_u_b_ref';
                'n_u_c_ref';
                'n_l_a_ref';
                'n_l_b_ref';
                'n_l_c_ref';
                };
    lasi_delayed = {
                'n_u_a_ref_dlyd';
                'n_u_b_ref_dlyd';
                'n_u_c_ref_dlyd';
                'n_l_a_ref_dlyd';
                'n_l_b_ref_dlyd';
                'n_l_c_ref_dlyd';
                };
    smsi_inputs = {
                'v_g_a';
                'v_g_b';
                'v_g_c';
                'v_d';
                };
    smsi_outputs = {
                'i_s_a';
                'i_s_b';
                'i_s_c';
                'i_d';
                };
    
    % _____________________________________________________________________
    %                                                  edit above this line
    
    DSA.vars.lasi_states   = lasi_states;
    DSA.vars.lasi_inputs   = lasi_inputs;
    DSA.vars.lasi_outputs  = lasi_outputs;
    DSA.vars.lasi_todelay  = lasi_todelay;
    DSA.vars.lasi_delayed  = lasi_delayed;
    DSA.vars.smsi_inputs   = smsi_inputs;
    DSA.vars.smsi_outputs  = smsi_outputs;
end