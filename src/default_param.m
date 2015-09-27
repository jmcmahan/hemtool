function dflt = default_param()

    % Type of function involved
    dflt.output_type = 'strain';
    dflt.function_type = 'forward';


    % Number of grid points for coercive and interaction fields
    dflt.Nc = 80;
    dflt.Ni = 80;


    % Remanence polarization
    dflt.P_r_p =  0.1773;
    dflt.P_r_m = -0.1773;
    dflt.P_r_n =  0.0;

    % Remanence strain
    dflt.eps_r_p = 4.638e-4;
    dflt.eps_r_m = 4.638e-4;
    dflt.eps_r_n = -2e-4;

    dflt.chi_p = 5e-8;
    dflt.chi_m = 5e-8;
    dflt.chi_n = 5e-8;

    dflt.d_p =  3.74e-10;
    dflt.d_m = -3.74e-10;
    dflt.d_n = 0.0;

    dflt.s_p = 1.88e-11;
    dflt.s_m = 1.88e-11;
    dflt.s_n = 1.88e-11;

    dflt.prestress = -10.6e6;

    % Distribution parameters
    dflt.ec_bar = 1.2e5;
    dflt.sigma_c = 2.554e-1;
    dflt.sigma_i = 2.5e5;

    % Weighted normal / lognormal distribution parameters
    dflt.sigma_i_scale = [1/4, 1/2, 1];
    dflt.sigma_c_scale = [1/4, 1/2, 1];
    dflt.mean_scale = [1];
    dflt.dist_i_weights = [1/10, 1/10, 1];
    % should have length(mean_scale)*length(sigma_c_scale) points for this
    dflt.dist_c_weights = [1/10, 1/10, 1];


    % Relaxation parameters
    dflt.gamma = 6e-4;
    dflt.tau = 7e-2;

    % Parameters specific to the 180 degree polarization model implementation
    dflt.beta = 3.7499999999999999e-03; 
    dflt.P_r = 0.1773;
    dflt.eta = 1.0308396819793936e7;

    % Input discretization parameters
    dflt.ecal_min = -2e6;
    dflt.ecal_max = 2e6;
    dflt.ecal_n = 257;
    dflt.delta_t = 1e-1;


    % Actuator parameters
    dflt.rod_area = 4.225e-5;
    dflt.rod_length = 1e-2;
    dflt.rod_stiffness = 2.7e6;

    % Quadrature setting
    dflt.quadrature = 'quad_4pt_gauss';


    % Distribution setting
    dflt.distribution = 'dist_lognormal';


    % Dipole evolution setting
    dflt.evolution = 'evol_90_relax';


    % HEM Evaluation setting
    dflt.evaluation = 'eval_90_strain';


    % Inverse function settings
    dflt.inverse = 'inverse_binary_search';
    dflt.search_range = [0, dflt.ecal_max];
    % This is just set to some low number for the strain. For polarization
    % you can set this to delta_E / (2*eta) to have the search exit when
    % the input has been found to within the exact tolerance.
    dflt.inverr_atol = 1e-6;
    % Scale the error by this value and add to result.
    dflt.correction_scale = 0;
    dflt.error_scale = 1;

end
