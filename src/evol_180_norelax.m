function hem = evol_180_norelax(hem, param)
    % Set up the update of the dipole fraction

    % Calculate the lookup tables used for fast evaluation of the dipole update

    % Field set up
    eim = hem.eim;
    ecm = hem.ecm;

    % Material parameters
    eta = param.eta;

    size_grid = size(eim); 
    m = size_grid(1);
    n = size_grid(2);

    % Save important parameters
    hem.eta = eta;
    hem.P_r = param.P_r;


    % Initialize dipole fraction variables
    hem.xp = ones(m, n);
    hem.xp(:,1:n/2) = 0;

    % Setup pointers to the update function
    hem.update = @(hem, E) evol_180_norelax_update(hem, E);
end


function [xp] = evol_180_norelax_update(hem, E)
    % This function is called (via hem.update) to find updated values for
    % the positve (xp) and negative (xm) dipole fractions.
    xp = hem.xp;

    eim = hem.eim;
    ecm = hem.ecm;

    % Here we converts 0->-1, 1->1, use the result to quickly update the
    % dipole fraction, and convert back. 
    xp = (sign(eim + E + ecm.*(2*xp - 1)) + 1) / 2;
end
