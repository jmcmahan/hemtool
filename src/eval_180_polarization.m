function hem = eval_180_polarization(hem, param)
    % Set up the functions for computing the polarization from the input 
    % and the updated dipole fraction.


    hem.integrate = @(hem, E) integrate_180_polarization(hem, E);
    hem.forward = @(hem, E) forward_180_polarization(hem, E);
    hem.forward_fix_t = @(hem, E) forward_180_polarization_fix_t(hem, E);

    hem.get_polarization = @(hem, E) forward_180_polarization(hem, E);
    hem.get_polarization_fix_t = @(hem, E) ...
                                    forward_180_polarization_fix_t(hem, E);
end

function [P, hem] = integrate_180_polarization(hem, E)
    % The evaluation of the integral and computation of the final output
    % occurs here.

    doquad = hem.doquad;

    eta = hem.eta;
    Pr = hem.P_r;

    xp = hem.xp;

    P = E/eta - Pr + 2*Pr*doquad(hem, xp);
end


function [P, hem] = forward_180_polarization(hem, E)
    % Computes the forward model for the input(s) E, updating the state
    % as it goes

    N = length(E);
    P = zeros(N, 1);
    for k = 1:N
        Ek = E(k);
        hem.xp = hem.update(hem, Ek);
        [P(k), hem] = hem.integrate(hem, Ek);
    end
end


function [P, hem] = forward_180_polarization_fix_t(hem, E)
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
        h.xp = hem.update(hem, Ek);
        P(k) = hem.integrate(h, Ek);
    end
end
