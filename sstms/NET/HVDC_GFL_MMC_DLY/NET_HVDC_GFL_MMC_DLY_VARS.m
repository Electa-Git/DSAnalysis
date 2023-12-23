function DSA = NET_HVDC_GFL_MMC_DLY_VARS(DSA)
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
    
    lasi_states = {
                % -- onshore MMC
                'i_s_1_al';
                'i_s_1_be';
                'i_c_1_a';
                'i_c_1_b';
                'i_c_1_c';
                'v_c_u_1_a';
                'v_c_u_1_b';
                'v_c_u_1_c';
                'v_c_l_1_a';
                'v_c_l_1_b';
                'v_c_l_1_c';
                'e_g_1_al_flt';
                'e_g_1_be_flt';
                'v_g_1_al_flt';
                'v_g_1_be_flt';
                'e_g_1_al_flt_quad';
                'e_g_1_be_flt_quad';
                'v_g_1_al_flt_quad';
                'v_g_1_be_flt_quad';
                'e_PLL_v_1';
                'theta_eps_v_1';
                'e_MV_1';
                'e_MQ_1';
                'e_MAC1_1_al';
                'e_MAC1_1_be';
                'e_MAC2_1_al';
                'e_MAC2_1_be';
                'e_MCC1_1_a';
                'e_MCC1_1_b';
                'e_MCC1_1_c';
                'e_MCC2_1_a';
                'e_MCC2_1_b';
                'e_MCC2_1_c';
                % -- offshore MMC
                'i_s_2_al';
                'i_s_2_be';
                'i_c_2_a';
                'i_c_2_b';
                'i_c_2_c';
                'v_c_u_2_a';
                'v_c_u_2_b';
                'v_c_u_2_c';
                'v_c_l_2_a';
                'v_c_l_2_b';
                'v_c_l_2_c';
                'e_g_2_al_flt';
                'e_g_2_be_flt';
                'v_g_2_al_flt';
                'v_g_2_be_flt';
                'e_g_2_al_flt_quad';
                'e_g_2_be_flt_quad';
                'v_g_2_al_flt_quad';
                'v_g_2_be_flt_quad';
                'e_PLL_v_2';
                'theta_eps_v_2';
                'e_MP_2';
                'e_MQ_2';
                'e_MAC1_2_al';
                'e_MAC1_2_be';
                'e_MAC2_2_al';
                'e_MAC2_2_be';
                'e_MCC1_2_a';
                'e_MCC1_2_b';
                'e_MCC1_2_c';
                'e_MCC2_2_a';
                'e_MCC2_2_b';
                'e_MCC2_2_c';
                % -- HVDC cables
                'v_d_u_1';
                'v_d_l_1';
                'i_z_u';
                'i_z_l';
                'v_d_u_2';
                'v_d_l_2';
                };
    lasi_inputs = {
                'v_g_1_a';
                'v_g_1_b';
                'v_g_1_c';
                'v_g_2_a';
                'v_g_2_b';
                'v_g_2_c';
                'v_d_1_ref';
                'q_g_1_ref';
                'p_g_2_ref';
                'q_g_2_ref';
                'cos_omega_1_t';
                'sin_omega_1_t';
                };
    lasi_outputs = {
                % -- HVDC cable
                'v_d_1';
                'v_d_del_1';
                'v_d_2';
                'v_d_del_2';
                % -- onshore MMC
                'i_s_1_a';
                'i_s_1_b';
                'i_s_1_c';
                'v_g_1_al';
                'v_g_1_be';
                'i_u_1_a';
                'i_u_1_b';
                'i_u_1_c';
                'i_l_1_a';
                'i_l_1_b';
                'i_l_1_c';
                'i_d_u_1';
                'i_d_l_1';
                'i_d_1';
                'i_d_del_1';
                'cos_theta_eps_1';
                'sin_theta_eps_1';
                'cos_theta_1';
                'sin_theta_1';
                'v_g_1_d_flt';
                'v_g_1_q_flt';
                'delta_theta_1';
                'delta_omega_1';
                'p_g_1';
                'q_g_1';
                'w_z_1_ref';
                'w_z_1';
                'i_s_1_d_ref'    
                'i_s_1_q_ref';
                'i_s_1_al_ref';
                'i_s_1_be_ref';
                'v_s_1_al_ref';
                'v_s_1_be_ref';
                'v_s_1_a_ref';
                'v_s_1_b_ref';
                'v_s_1_c_ref';
                'i_c_1_ref';
                'v_c_1_a_ref';
                'v_c_1_b_ref';
                'v_c_1_c_ref';
                'n_u_1_a_ref';
                'n_u_1_b_ref';
                'n_u_1_c_ref';
                'n_l_1_a_ref';
                'n_l_1_b_ref';
                'n_l_1_c_ref';
                'n_u_1_a';
                'n_u_1_b';
                'n_u_1_c';
                'n_l_1_a';
                'n_l_1_b';
                'n_l_1_c';
                'v_u_1_a';
                'v_u_1_b';
                'v_u_1_c';
                'v_l_1_a';
                'v_l_1_b';
                'v_l_1_c';
                'v_s_1_a';
                'v_s_1_b';
                'v_s_1_c';
                'v_c_1_a';
                'v_c_1_b';
                'v_c_1_c';
                'v_s_1_al';
                'v_s_1_be';
                % -- offshore MMC
                'i_s_2_a';
                'i_s_2_b';
                'i_s_2_c';
                'v_g_2_al';
                'v_g_2_be';
                'i_u_2_a';
                'i_u_2_b';
                'i_u_2_c';
                'i_l_2_a';
                'i_l_2_b';
                'i_l_2_c';
                'i_d_u_2';
                'i_d_l_2';
                'i_d_2';
                'i_d_del_2';
                'cos_theta_eps_2';
                'sin_theta_eps_2';
                'cos_theta_2';
                'sin_theta_2';
                'v_g_2_d_flt';
                'v_g_2_q_flt';
                'delta_theta_2';
                'delta_omega_2';
                'p_g_2';
                'q_g_2';
                'i_s_2_d_ref';
                'i_s_2_q_ref';
                'i_s_2_al_ref';
                'i_s_2_be_ref';
                'v_s_2_al_ref';
                'v_s_2_be_ref';
                'v_s_2_a_ref';
                'v_s_2_b_ref';
                'v_s_2_c_ref';
                'i_c_2_ref';
                'v_c_2_a_ref';
                'v_c_2_b_ref';
                'v_c_2_c_ref';
                'n_u_2_a_ref';
                'n_u_2_b_ref';
                'n_u_2_c_ref';
                'n_l_2_a_ref';
                'n_l_2_b_ref';
                'n_l_2_c_ref';
                'n_u_2_a';
                'n_u_2_b';
                'n_u_2_c';
                'n_l_2_a';
                'n_l_2_b';
                'n_l_2_c';
                'v_u_2_a';
                'v_u_2_b';
                'v_u_2_c';
                'v_l_2_a';
                'v_l_2_b';
                'v_l_2_c';
                'v_s_2_a';
                'v_s_2_b';
                'v_s_2_c';
                'v_c_2_a';
                'v_c_2_b';
                'v_c_2_c';
                'v_s_2_al';
                'v_s_2_be';
                };
    lasi_todelay = {
                'n_u_1_a_ref';
                'n_u_1_b_ref';
                'n_u_1_c_ref';
                'n_l_1_a_ref';
                'n_l_1_b_ref';
                'n_l_1_c_ref';
                'n_u_2_a_ref';
                'n_u_2_b_ref';
                'n_u_2_c_ref';
                'n_l_2_a_ref';
                'n_l_2_b_ref';
                'n_l_2_c_ref';
                };
    lasi_delayed = {
                'n_u_1_a_ref_dlyd';
                'n_u_1_b_ref_dlyd';
                'n_u_1_c_ref_dlyd';
                'n_l_1_a_ref_dlyd';
                'n_l_1_b_ref_dlyd';
                'n_l_1_c_ref_dlyd';
                'n_u_2_a_ref_dlyd';
                'n_u_2_b_ref_dlyd';
                'n_u_2_c_ref_dlyd';
                'n_l_2_a_ref_dlyd';
                'n_l_2_b_ref_dlyd';
                'n_l_2_c_ref_dlyd';
                };
    smsi_inputs = {
                'v_g_1_a';
                'v_g_1_b';
                'v_g_1_c';
                'v_g_2_a';
                'v_g_2_b';
                'v_g_2_c';
                };
    smsi_outputs = {
                'i_s_1_a';
                'i_s_1_b';
                'i_s_1_c';
                'i_s_2_a';
                'i_s_2_b';
                'i_s_2_c';
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