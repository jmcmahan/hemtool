function hem = eval_90_polarization(hem, param)
    % Set up the functions for computing the polarization from the input 
    % and the updated dipole fraction. This particular form of 
    % polarization uses the strain computation, so it's initialized
    % first

    hem = eval_90_strain(hem, param);
    % We save this since we use it when updating the polarization
    hem.integrate_90_strain = hem.integrate;

    hem.integrate = @(hem, E) integrate_90_polarization(hem, E);
    hem.forward = @(hem, E) forward_90_polarization(hem, E);
    hem.forward_fix_t = @(hem, E) forward_90_polarization_fix_t(hem, E);

    hem.get_polarization = @(hem, E) forward_90_polarization(hem, E);
    hem.get_polarization_fix_t = @(hem, E) ...
                                    forward_90_polarization_fix_t(hem, E);
end

function [P, hem] = integrate_90_polarization(hem, E)
    % The evaluation of the integral and computation of the final output
    % occurs here.

    disti = hem.disti;
    distc = hem.distc;

    xp = hem.xp;
    xm = hem.xm;

    i = hem.get_input_index(hem, E);

    doquad = hem.doquad;

    P_m_p = hem.P_m_p(:,:,i);
    P_m_m = hem.P_m_m(:,:,i);
    P_m_n = hem.P_m_n(:,:,i);

    P_bar = xp.*(P_m_p - P_m_n) + xm.*(P_m_m - P_m_n) + P_m_n;

    P = doquad(hem, P_bar);
end


function [P, hem] = forward_90_polarization(hem, E)
    % Computes the forward model for the input(s) E, updating the state
    % as it goes

    N = length(E);
    P = zeros(N, 1);
    for k = 1:N
        Ek = E(k);
        [hem.xp, hem.xm] = hem.update(hem, Ek);
        [P(k), hem] = hem.integrate(hem, Ek);
    end
end


function [P, hem] = forward_90_polarization_fix_t(hem, E)
    % This version of the function does not update the dipole fraction as
    % it evaluates the strain at each E, so this it gives the strain for
    % the various possibilities of E for a fixed point in time.

    N = length(E);
    P = zeros(N, 1);
    h = hem;
    for k = 1:N
        Ek = E(k);
        % Here we are simply discarding the returned updated structure, which
        % has the effect of fixing the state at the current time.
        [h.xp, h.xm] = hem.update(hem, Ek);
        P(k) = hem.integrate(h, Ek);
    end
end
