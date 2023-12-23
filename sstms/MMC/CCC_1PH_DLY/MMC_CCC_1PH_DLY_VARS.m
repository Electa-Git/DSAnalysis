function DSA = MMC_CCC_1PH_DLY_VARS(DSA)
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
                'v_g';
                'v_d';
                'i_s_ref';
                };
    lasi_states = {
                'i_s';
                'i_c';
                'v_c_u';
                'v_c_l';
                'i_d_flt';
                'eta_AC1';
                'eta_AC2';
                'eta_CC1';
                'eta_CC2';
                };
    lasi_outputs = {
                'i_s';
                'i_u';
                'i_l';
                'i_c_ref';
                'v_s_ref';
                'v_c_ref';
                'v_s';
                'v_c';
                'p_g';
                'i_d';
                'n_u';
                'n_l';
                'v_u';
                'v_l';
                };
    lasi_todelay = {
                'n_u_ref';
                'n_l_ref';
                };
    lasi_delayed = {
                'n_u';
                'n_l';
                };
    smsi_inputs = {
                'v_g';
                'v_d';
                };
    smsi_outputs = {
                'i_s';
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