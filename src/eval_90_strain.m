function hem = eval_90_strain(hem, param)
    % Set up the functions for computing the strain from the input and the
    % updated dipole fraction.

    A = param.rod_area;
    L = param.rod_length;
    ks = param.rod_stiffness;

    hem.rod_coeff = L * ks / A;

    hem.integrate = @(hem, E) integrate_90_strain(hem, E);
    hem.forward = @(hem, E) forward_90_strain(hem, E);
    hem.forward_fix_t = @(hem, E) forward_90_strain_fix_t(hem, E);

    hem.get_strain = @(hem, E) forward_90_strain(hem, E);
    hem.get_strain_fix_t = @(hem, E) forward_90_strain_fix_t(hem, E);
end

function [eps, hem] = integrate_90_strain(hem, E)
    % The evaluation of the integral and computation of the final output
    % occurs here.

    doquad = hem.doquad;

    sp = hem.s_p;
    sm = hem.s_m;
    sn = hem.s_n;
    s = hem.prestress;
    dp = hem.d_p;
    epsrp = hem.eps_r_p;
    epsrn = hem.eps_r_n;

    xp = hem.xp;
    xm = hem.xm;


    %i = hem.get_input_index(hem, E);

    %eps_m_p = hem.eps_m_p(:,:,i);
    %eps_m_m = hem.eps_m_m(:,:,i);
    %eps_m_n = hem.eps_m_n(:,:,i);

    %eps_bar = xp.*(eps_m_p - eps_m_n) + xm.*(eps_m_m - eps_m_n) + eps_m_n;
    eps_bar = E*dp*(xp - xm) + (epsrp - epsrn)*(xp - xm);

    eps = s*sp + epsrn + doquad(hem, eps_bar);
end


function [eps, hem] = forward_90_strain(hem, E)
    % Computes the forward model for the input(s) E, updating the state
    % as it goes

    N = length(E);
    eps = zeros(N, 1);
    for k = 1:N
        Ek = E(k);
        [hem.xp, hem.xm] = hem.update(hem, Ek);
        [eps(k), hem] = hem.integrate(hem, Ek);
    end
end


function [eps, hem] = forward_90_strain_fix_t(hem, E)
    % This version of the function does not update the dipole fraction as
    % it evaluates the strain at each E, so this it gives the strain for
    % the various possibilities of E for a fixed point in time.

    N = length(E);
    eps = zeros(N, 1);
    h = hem;
    for k = 1:N
        Ek = E(k);
        % Here we are simply discarding the returned updated structure, which
        % has the effect of fixing the state at the current time.
        [h.xp, h.xm] = hem.update(hem, Ek);
        eps(k) = hem.integrate(h, Ek);
    end
end
