function DSA = NET_GFM_MMC_OWF_3TLC_DLY_VARS(DSA)
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
    % /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
    % The order in which the variables are listed is important:
    %     - The state variables must be in the same order as in the dxdt
    %      vector in the EQNS file.
    %     - The delayed and todelay lists must have matching order.
    % /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
    % 
    % The only difference between lasi and smsi systems are the
    % lists of inputs and outputs. The states, todelay and delayed
    % variables remain the same.
    % _____________________________________________________________________
    %                                                  edit below this line
    
    lasi_inputs = {
                'v_d_2';
                'v_g_2_a_ref';
                'v_g_2_b_ref';
                'v_g_2_c_ref';
                'p_f_1_ref';
                'q_f_1_ref';
                'p_f_2_ref';
                'q_f_2_ref';
                'p_f_3_ref';
                'q_f_3_ref';
                'v_d_w_1';
                'v_d_w_2';
                'v_d_w_3';
                'cos_omega_1_t';
                'sin_omega_1_t';
                };
    lasi_states = {
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
                'e_MAV1_al';
                'e_MAV1_be';
                'e_MAV2_al';
                'e_MAV2_be';
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
                % -- offshore PCC
                'v_g_2_a';
                'v_g_2_b';
                'v_g_2_c';
                % -- TCL 1 and collection cable 1
                'i_f_1_al';
                'i_f_1_be';
                'v_f_1_al';
                'v_f_1_be';
                'i_t_1_al';
                'i_t_1_be';
                'v_x_1_al';
                'v_x_1_be';
                'i_x_1_al';
                'i_x_1_be';
                'e_f_1_al_flt';
                'e_f_1_be_flt';
                'v_f_1_al_flt';
                'v_f_1_be_flt';
                'e_f_1_al_flt_quad';
                'e_f_1_be_flt_quad';
                'v_f_1_al_flt_quad';
                'v_f_1_be_flt_quad';
                'e_PLL_w_1';
                'theta_eps_w_1';
                'e_TPQ_1_d';
                'e_TPQ_1_q';
                'e_TAC1_1_al';
                'e_TAC1_1_be';
                'e_TAC2_1_al';
                'e_TAC2_1_be';
                % -- TCL 2 and collection cable 2
                'i_f_2_al';
                'i_f_2_be';
                'v_f_2_al';
                'v_f_2_be';
                'i_t_2_al';
                'i_t_2_be';
                'v_x_2_al';
                'v_x_2_be';
                'i_x_2_al';
                'i_x_2_be';
                'e_f_2_al_flt';
                'e_f_2_be_flt';
                'v_f_2_al_flt';
                'v_f_2_be_flt';
                'e_f_2_al_flt_quad';
                'e_f_2_be_flt_quad';
                'v_f_2_al_flt_quad';
                'v_f_2_be_flt_quad';
                'e_PLL_w_2';
                'theta_eps_w_2';
                'e_TPQ_2_d';
                'e_TPQ_2_q';
                'e_TAC1_2_al';
                'e_TAC1_2_be';
                'e_TAC2_2_al';
                'e_TAC2_2_be';
                % -- TCL 3 and collection cable 3
                'i_f_3_al';
                'i_f_3_be';
                'v_f_3_al';
                'v_f_3_be';
                'i_t_3_al';
                'i_t_3_be';
                'v_x_3_al';
                'v_x_3_be';
                'i_x_3_al';
                'i_x_3_be';
                'e_f_3_al_flt';
                'e_f_3_be_flt';
                'v_f_3_al_flt';
                'v_f_3_be_flt';
                'e_f_3_al_flt_quad';
                'e_f_3_be_flt_quad';
                'v_f_3_al_flt_quad';
                'v_f_3_be_flt_quad';
                'e_PLL_w_3';
                'theta_eps_w_3';
                'e_TPQ_3_d';
                'e_TPQ_3_q';
                'e_TAC1_3_al';
                'e_TAC1_3_be';
                'e_TAC2_3_al';
                'e_TAC2_3_be';
                };
    lasi_outputs = {
                % -- offshore MMC
                'i_s_2_a';
                'i_s_2_b';
                'i_s_2_c';
                'i_s_2_a_ref';
                'i_s_2_b_ref';
                'i_s_2_c_ref';
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
                'p_g_2';
                'q_g_2';
                'v_s_2_a_ref';
                'v_s_2_b_ref';
                'v_s_2_c_ref';
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
                % -- TCL 1 and collection cable 1
                'i_f_1_a';
                'i_f_1_b';
                'i_f_1_c';
                'v_t_1_al';
                'v_t_1_be';
                'i_y_1_al';
                'i_y_1_be';
                'i_x_1_a';
                'i_x_1_b';
                'i_x_1_c';
                'v_f_1_al_flt_ps';
                'v_f_1_be_flt_ps';
                'v_f_1_al_flt_ns';
                'v_f_1_be_flt_ns';
                'cos_theta_w_1';
                'sin_theta_w_1';
                'v_f_1_d_flt_ps';
                'v_f_1_q_flt_ps';
                'p_f_1';
                'q_f_1';
                'i_f_1_d_ref';
                'i_f_1_q_ref';
                'i_f_1_al_ref';
                'i_f_1_be_ref';
                'v_w_1_al_ref';
                'v_w_1_be_ref';
                'm_w_1_al_ref';
                'm_w_1_be_ref';
                'm_w_1_al';
                'm_w_1_be';
                'v_w_1_al';
                'v_w_1_be';
                % -- TCL 2 and collection cable 2
                'i_f_2_a';
                'i_f_2_b';
                'i_f_2_c';
                'v_t_2_al';
                'v_t_2_be';
                'i_y_2_al';
                'i_y_2_be';
                'i_x_2_a';
                'i_x_2_b';
                'i_x_2_c';
                'v_f_2_al_flt_ps';
                'v_f_2_be_flt_ps';
                'v_f_2_al_flt_ns';
                'v_f_2_be_flt_ns';
                'cos_theta_w_2';
                'sin_theta_w_2';
                'v_f_2_d_flt_ps';
                'v_f_2_q_flt_ps';
                'p_f_2';
                'q_f_2';
                'i_f_2_d_ref';
                'i_f_2_q_ref';
                'i_f_2_al_ref';
                'i_f_2_be_ref';
                'v_w_2_al_ref';
                'v_w_2_be_ref';
                'm_w_2_al_ref';
                'm_w_2_be_ref';
                'm_w_2_al';
                'm_w_2_be';
                'v_w_2_al';
                'v_w_2_be';
                % -- TCL 3 and collection cable 3
                'i_f_3_a';
                'i_f_3_b';
                'i_f_3_c';
                'v_t_3_al';
                'v_t_3_be';
                'i_y_3_al';
                'i_y_3_be';
                'i_x_3_a';
                'i_x_3_b';
                'i_x_3_c';
                'v_f_3_al_flt_ps';
                'v_f_3_be_flt_ps';
                'v_f_3_al_flt_ns';
                'v_f_3_be_flt_ns';
                'cos_theta_w_3';
                'sin_theta_w_3';
                'v_f_3_d_flt_ps';
                'v_f_3_q_flt_ps';
                'p_f_3';
                'q_f_3';
                'i_f_3_d_ref';
                'i_f_3_q_ref';
                'i_f_3_al_ref';
                'i_f_3_be_ref';
                'v_w_3_al_ref';
                'v_w_3_be_ref';
                'm_w_3_al_ref';
                'm_w_3_be_ref';
                'm_w_3_al';
                'm_w_3_be';
                'v_w_3_al';
                'v_w_3_be';
                % -- offshore PCC
                'i_g_2_a';
                'i_g_2_b';
                'i_g_2_c';
                'i_x_sum_a';
                'i_x_sum_b';
                'i_x_sum_c';
                };
    lasi_todelay = {
                'n_u_2_a_ref';
                'n_u_2_b_ref';
                'n_u_2_c_ref';
                'n_l_2_a_ref';
                'n_l_2_b_ref';
                'n_l_2_c_ref';
                'm_w_1_al_ref';
                'm_w_1_be_ref';
                'm_w_2_al_ref';
                'm_w_2_be_ref';
                'm_w_3_al_ref';
                'm_w_3_be_ref';
                };
    lasi_delayed = {
                'n_u_2_a_ref_dlyd';
                'n_u_2_b_ref_dlyd';
                'n_u_2_c_ref_dlyd';
                'n_l_2_a_ref_dlyd';
                'n_l_2_b_ref_dlyd';
                'n_l_2_c_ref_dlyd';
                'm_w_1_al_ref_dlyd';
                'm_w_1_be_ref_dlyd';
                'm_w_2_al_ref_dlyd';
                'm_w_2_be_ref_dlyd';
                'm_w_3_al_ref_dlyd';
                'm_w_3_be_ref_dlyd';
                };
    smsi_inputs = {
                'v_d_2';
                'v_d_w_1';
                };
    smsi_outputs = {
                'i_d_2';
                'i_d_w_1';
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