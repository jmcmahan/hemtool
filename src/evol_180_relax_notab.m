function hem = evol_180_relax_notab(hem, param)
    % Set up the update of the dipole fraction

    % Calculate the lookup tables used for fast evaluation of the dipole update

    % Field set up

    ecal_min = param.ecal_min;
    ecal_max = param.ecal_max;
    N = param.ecal_n;
    delta_E = (ecal_max - ecal_min) / (N - 1);
    ecal = ecal_min:delta_E:ecal_max;

    eim = hem.eim;

    eta = param.eta;

    tau = param.tau;
    beta = param.beta;
    dt = param.delta_t;
    gamma1 = sqrt(2 / pi) / (beta * tau);
    gamma2 = 1 / (sqrt(2) * beta * eta);



    % Save important parameters
    hem.eta = eta;
    hem.P_r = param.P_r;
    hem.tau = tau;
    hem.gamma1 = gamma1;
    hem.gamma2 = gamma2;
    hem.beta = beta;

    % Save computed values

    hem.ecal = ecal;
    hem.ecal_min = ecal_min;
    hem.ecal_max = ecal_max;
    hem.ecal_n = N;
    hem.delta_E = delta_E;


    hem.param = param;
    % Initialize dipole fraction variables
    size_grid = size(eim);
    m = size_grid(1);
    n = size_grid(2);
    hem.xp = ones(m,n);
    hem.xp(:,1:n/2) = 0;

    % Setup pointers to the update function
    hem.update = @(hem, E) evol_180_relax_update(hem, E);
end




function [xp] = evol_180_relax_update(hem, E)
    % This function is called (via hem.update) to find updated values for
    % the positve (xp) and negative (xm) dipole fractions.
    param = hem.param;
    xp = hem.xp;

    eim = hem.eim;
    ecm = hem.ecm;

    % Material parameters
    eta = param.eta;

    tau = param.tau;
    beta = param.beta;
    dt = param.delta_t;
    gamma1 = sqrt(2 / pi) / (beta * tau);
    gamma2 = 1 / (sqrt(2) * beta * eta);

    E = E + eim;

    % C is the value which acts on X = [x_+ x_-]' and D is the value 
    % added to the system when applying implicit Euler 
    % (i.e., X_next = C*X_cur + D)

    p_m_p = gamma1 ./ erfcx(-gamma2*(ecm - E));
    p_p_m = gamma1 ./ erfcx(-gamma2*(ecm + E));

    C = 1 ./ (1 + dt*(p_p_m + p_m_p));
    D = dt * p_m_p .* C;

    %C = hem.C_eul(:,:,i);
    %D = hem.D_eul(:,:,i);

    xp = C.*xp +  D;
end
