function hem = inverse_binary_search(hem, param)
    % Set up the function for performing the approximate hysteresis inversion
    % via a binary search.

    N = param.ecal_n;
    ecal_max = param.ecal_max;
    ecal_min = param.ecal_min;
    hem.delta_E = (ecal_max - ecal_min) / (N - 1);
    hem.search_range = [ecal_min ecal_max];
    hem.search_range = param.search_range;
    hem.inverr_atol = param.inverr_atol;
    hem.correction_scale = param.correction_scale;
    hem.error_scale = param.error_scale;
    hem.invert = @(out_d, hem) do_binary_search(out_d, hem);
end


function [E, hem] = do_binary_search(hem, out_d)

    % 'delta_E' is only naturally defined when using the discretized input, so
    % if using with some other form of the forward model (e.g., neglecting
    % relaxation), this needs to be defined when initializing 'hem'.
    delta_E = hem.delta_E; 
    emin = hem.search_range(1);
    emax = hem.search_range(2);
    atol = hem.inverr_atol;
    dEmin = delta_E / 4;
    dEmin = delta_E / 16;
    cs = hem.correction_scale;
    es = hem.error_scale;

    N = length(out_d);
    E = zeros(N,1);

    for k = 1:N
        Ek = emin + (emax - emin) / 2;
        dE = (emax - emin) / 4;
        % Note we save the output as hc, rather than updating hem
        [out_k, hc] = hem.forward(hem, Ek);
        err_k = es * (out_d(k) - out_k);
        while (abs(err_k) >= atol) & (dE >= dEmin)
            Ek = Ek + sign(err_k)*dE;
            [out_k, hc] = hem.forward(hem, Ek);
            err_k = es * (out_d(k) - out_k);
            dE = dE / 2;
        end
        E(k) = Ek + cs*err_k;
        % Now that the correct value has been found, we updated hem, which
        % updates the state for the next step.
        hem = hc;
    end
end
