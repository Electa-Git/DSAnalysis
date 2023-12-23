function HDi = RLC_Ch3_HD(DSA, HDi)
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
    HDi.f_b = 1; % [Hz]
    
    % extracting parameters:
    V_s = param.V_s;
    
    % --- Input variables:
    
    % defining a constant component (harmonic 0) as well as 
    % fundamental-frequency component (harmonics +/-1)
    HDi.v_s(h2i(0))      = V_s;
    HDi.v_s(h2i([-1 1])) = DSA_sine(0.5*V_s, 0);
    
    % --- State variables:
    % (just initial guesses)
    
    % if not specified here, automatically initialised to zero by the code.
    
    % --- Delayed variables:
    % (just initial guesses)
    
    % none
end