function hem = evol_90_relax(hem, param)
    % Set up the update of the dipole fraction

    % Calculate the lookup tables used for fast evaluation of the dipole update

    % Field set up

    ecal_min = param.ecal_min;
    ecal_max = param.ecal_max;
    N = param.ecal_n;
    delta_E = (ecal_max - ecal_min) / (N - 1);
    ecal = ecal_min:delta_E:ecal_max;

    eim = hem.eim;
    ecm = hem.ecm;

    % Material parameters
    P_r_p = param.P_r_p;
    P_r_m = param.P_r_m;
    P_r_n = param.P_r_n;

    eps_r_p = param.eps_r_p;
    eps_r_m = param.eps_r_m;
    eps_r_n = param.eps_r_n;

    chi_p = param.chi_p;
    chi_m = param.chi_m;
    chi_n = param.chi_n;

    d_p = param.d_p;
    d_m = param.d_m;
    d_n = param.d_n;

    s_p = param.s_p;
    s_m = param.s_m;
    s_n = param.s_n;

    s = param.prestress;

    gamma = param.gamma;
    tau = param.tau;
    dt = param.delta_t;

    %dg0 = param.ec_bar / 4;
    dg0 = ecm / 4;

    size_grid = size(eim); 
    m = size_grid(1);
    n = size_grid(2);
    p = N;
    C11_eul = zeros(m, n, p);
    C12_eul = zeros(m, n, p);
    C21_eul = zeros(m, n, p);
    C22_eul = zeros(m, n, p);
    D1_eul = zeros(m, n, p);
    D2_eul = zeros(m, n, p);
    eps_m_p = zeros(m, n, p);
    eps_m_m = zeros(m, n, p);
    eps_m_n = zeros(m, n, p);
    P_m_p = zeros(m, n, p);
    P_m_m = zeros(m, n, p);
    P_m_n = zeros(m, n, p);

    for i = 1:N
        E = ecal(i) + eim;

        eps_m_p_i = eps_r_p + d_p*E + s_p*s;
        eps_m_m_i = eps_r_m + d_m*E + s_m*s;
        eps_m_n_i = eps_r_n + d_n*E + s_n*s;

        P_m_p_i = P_r_p + chi_p*E + d_p*s;
        P_m_m_i = P_r_m + chi_m*E + d_m*s;
        P_m_n_i = P_r_n + chi_n*E + d_n*s;

        G_m_p = -E.*(P_m_p_i + P_r_p) / 2 - s.*(eps_m_p_i + eps_r_p) / 2;
        G_m_m = -E.*(P_m_m_i + P_r_m) / 2 - s.*(eps_m_m_i + eps_r_m) / 2;
        G_m_n = -E.*(P_m_n_i + P_r_n) / 2 - s.*(eps_m_n_i + eps_r_n) / 2;

        f_p_n = G_m_p - G_m_n;
        f_n_m = G_m_n - G_m_m;
        f_m_n = -f_n_m;
        f_n_p = -f_p_n;

        dg_p_n = dg0.*(1 - f_p_n./ecm).^2 .* (f_p_n < ecm);
        dg_n_p = dg0.*(1 - f_n_p./ecm).^2 .* (f_n_p < ecm);
        dg_m_n = dg0.*(1 - f_m_n./ecm).^2 .* (f_m_n < ecm);
        dg_n_m = dg0.*(1 - f_n_m./ecm).^2 .* (f_n_m < ecm);

        p_p_n = exp(-dg_p_n * gamma) / tau;
        p_n_p = exp(-dg_n_p * gamma) / tau;
        p_m_n = exp(-dg_m_n * gamma) / tau;
        p_n_m = exp(-dg_n_m * gamma) / tau;

        % Numbers here corresponding to row / column indices of a 2x2 system
        a11 = 1 + dt*(p_p_n + p_n_p);
        a12 = dt * p_n_p;
        a22 = 1 + dt*(p_m_n + p_n_m); 
        a21 = dt * p_n_m;
        detA = a11.*a22 - a12.*a21;

        % C is the matrix which acts on X = [x_+ x_-]' and D is the vector 
        % added to the system when applying implicit Euler 
        % (i.e., X_next = C*X_cur + D)

        C11_eul(:,:,i) =  a22 ./ detA; 
        C12_eul(:,:,i) = -a12 ./ detA;
        C21_eul(:,:,i) = -a21 ./ detA; 
        C22_eul(:,:,i) =  a11 ./ detA;

        D1_eul(:,:,i) = ( a22.*a12 - a12.*a21) ./ detA;
        D2_eul(:,:,i) = (-a21.*a12 + a11.*a21) ./ detA; 

        eps_m_p(:,:,i) = eps_m_p_i;
        eps_m_m(:,:,i) = eps_m_m_i;
        eps_m_n(:,:,i) = eps_m_n_i;

        P_m_p(:,:,i) = P_m_p_i;
        P_m_m(:,:,i) = P_m_m_i;
        P_m_n(:,:,i) = P_m_n_i;

        % If you want to save the transition rates, can do so here
        %p_pn(:,:,i) = p_p_n;
        %p_np(:,:,i) = p_n_p;
        %p_mn(:,:,i) = p_m_n;
        %p_nm(:,:,i) = p_n_m;

    end

    % Save important parameters
    hem.prestress = s;
    hem.s_p = s_p;
    hem.s_m = s_m;
    hem.s_n = s_n;
    hem.d_p = d_p;
    hem.d_m = d_m;
    hem.d_n = d_n;
    hem.eps_r_p = eps_r_p;
    hem.eps_r_n = eps_r_n;

    % Save computed values

    hem.ecal = ecal;
    hem.ecal_min = ecal_min;
    hem.ecal_max = ecal_max;
    hem.ecal_n = N;
    hem.delta_E = delta_E;

    % This arrangement cuts down on arithmetic operations when evaluating
    hem.eps_m_p = eps_m_p;
    hem.eps_m_m = eps_m_m;
    hem.eps_m_n = eps_m_n;

    hem.P_m_p = P_m_p;
    hem.P_m_m = P_m_m;
    hem.P_m_n = P_m_n;

    hem.C11_eul = C11_eul;
    hem.C12_eul = C12_eul;
    hem.C21_eul = C21_eul;
    hem.C22_eul = C22_eul;
    hem.D1_eul = D1_eul;
    hem.D2_eul = D2_eul;

    %hem.p_pn = p_pn;
    %hem.p_np = p_np;
    %hem.p_mn = p_mn;
    %hem.p_nm = p_nm;

    % Initialize dipole fraction variables
    hem.xp = ones(m, n);
    hem.xm = zeros(m, n);
    hem.xp(:,1:n/2) = 0;
    hem.xm(:,1:n/2) = 1;

    % Setup pointers to the update function
    hem.get_input_index = @(hem, E) get_input_index(hem, E);
    hem.update = @(hem, E) evol_90_relax_update(hem, E);

end


function i = get_input_index(hem, E)
    % Round the input to the nearest value on the discrete input grid and
    % return the corresponding index.

    Emin = hem.ecal_min;
    dE = hem.delta_E;
    N = hem.ecal_n;

    i = round((E - Emin) / dE) + 1;
    i = max(1, i);
    i = min(N, i);
end


function [xp, xm] = evol_90_relax_update(hem, E)
    % This function is called (via hem.update) to find updated values for
    % the positve (xp) and negative (xm) dipole fractions.
    xp = hem.xp;
    xm = hem.xm;

    i = hem.get_input_index(hem, E);

    C11 = hem.C11_eul(:,:,i);
    C12 = hem.C12_eul(:,:,i);
    C21 = hem.C21_eul(:,:,i);
    C22 = hem.C22_eul(:,:,i);
    D1 = hem.D1_eul(:,:,i);
    D2 = hem.D2_eul(:,:,i);

    xp = C11.*xp + C12.*xm + D1;
    xm = C21.*xp + C22.*xm + D2;
end
