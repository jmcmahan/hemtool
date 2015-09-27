function hem = dist_lognormal(hem, param)
    % Computes the normal / lognormal distribution and merges with
    % the quadrature weights for easy computation.

    % These are vectors giving the coefficients to multiply the
    % paramaters, sigma_i, sigma_c, and sigma_i by when using a 
    % weighted sum of distributions. So e.g.
    % sigma_i_scale = [1, 1/2, 1/4]
    % would generate the densities for 
    % 1*sigma_i, 1/2*sigma_i, 1/4*sigma_i.
    % The values in 'dist_i_weights' give the weights these densities
    % are multiplied by before summing, so in this example,
    % dist_i_weights needs to be a vector of dimension 3.
    %
    % It is similar for sigma_c, except it has the additional parameter
    % 'mean_scale', so there will be 
    % length(sigma_c_scale) * length(mean_scale)
    % coefficients for the coercive distribution, which means there
    % will be this many values in "dist_c_weights" as well.
    sis = param.sigma_i_scale;
    scs = param.sigma_c_scale;
    ms = param.mean_scale;

    % Required: length(diw) = length(sis)
    diw = param.dist_i_weights;
    % Required: length(diw) = length(sis)*length(ms)
    dcw = param.dist_c_weights;

    si = param.sigma_i;
    sc = param.sigma_c;
    ec_bar = param.ec_bar;

    ei = hem.ei_gpts;
    ec = hem.ec_gpts;

    % These quadrature weights will be combined with the distributions
    % computed for easy evaluation of the integral.

    wi = hem.wi;
    wc = hem.wc;

    Ni = hem.Ni;
    Nc = hem.Nc;

    quadi = [];
    quadc = [];

    % Use these if you want to look at the distributions for whatever
    % reason
    nui = [];
    nuc = [];


    % Build normal interaction distribution combined with quadrature weights

    for i = 1:length(sis)
        s = si * sis(i);
        nu_i_pos = exp(-(ei(Ni/2 + 1: Ni) / s).^2 / 2);
        %nu_i_neg = rot90(rot90(nu_i_pos));
        nu_i_neg = flipud(nu_i_pos);
        nu_i = [nu_i_neg; nu_i_pos] / (sqrt(2*pi) * s);

        quadi = [quadi, nu_i .* wi];
        nui = [nui, nu_i];
    end

    % Build lognormal coercive distribution combined with quadrature weights

    for i = 1:length(scs)
        s = sc * scs(i);
        for j = 1:length(ms)
            mu = log(ms(j) * ec_bar);
            nu_c = exp(-( (log(ec) - mu).^2 / (2*s^2)));
            nu_c = nu_c ./ (ec * sqrt(2*pi) * s);
            quadc = [quadc, nu_c .* wc];
            nuc = [nuc, nu_c];
        end
    end


    % Note - for some applications, such as adaptive estimation of the
    % distribution, it may be useful to change the weights for the 
    % distribution. In those cases, rather than doing the following, you
    % would save "quadi" and "quadc" and do this multiplication on the
    % fly.

    % In evaluating the integral, distc multiplies on the left, 
    % disti on the right (no need to transpose since it's done here).
    hem.disti = quadi * diw' / sum(diw);
    hem.distc = dcw * quadc' / sum(dcw);

    % We'll go ahead and ensure everything evaluates to 1.
    intconst = sqrt(hem.distc * ones(Nc,Ni) * hem.disti);
    hem.disti = hem.disti / intconst;
    hem.distc = hem.distc / intconst;

    % Save these in case you want to look at the distributions without
    % the quadrature weights
    hem.nui = nui * diw';
    hem.nuc = nuc * dcw';

end
