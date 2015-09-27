function hem = quad_4pt_gauss(hem, param)
    % Calculates the grid points and quadrature weights for a 4-point
    % Gaussian quadrature rule. 

    % Nc and Ni should be divisible by 4.
    Nc = param.Nc;
    Ni = param.Ni;

    ec_bar = param.ec_bar;
    sigma_c = param.sigma_c;
    sigma_i = param.sigma_i;
        
    % Set integration limits
    if ~isfield(param, 'ec_int')
        ec_int = ec_bar*exp(4*sigma_c);
    else
        ec_int = param.ec_int;
    end

    if ~isfield(param, 'ei_int')
        ei_int = 4*sigma_i;
    else
        ei_int = param.ei_int;
    end

    % Setup coercive field grid points
    hc = 4 * ec_int / Nc;
    wc = gauss_weights(Nc/4, hc);
    ec_gpts = gauss_points(Nc/4, hc, 0);

    % Setup interaction field grid points
    hi = 8 * ei_int / Ni;
    wi = gauss_weights(Ni/4, hi);
    ei_gpts = gauss_points(Ni/4, hi, -ei_int);

    % Save results into structure
    hem.Nc = Nc;
    hem.ec_gpts = ec_gpts;
    hem.wc = wc;
    % Matrix of coercive grid points
    hem.ecm = ec_gpts*ones(1, Ni);

    hem.Ni = Ni;
    hem.ei_gpts = ei_gpts;
    hem.wi = wi;
    % Matrix of interaction grid points
    hem.eim = ones(Nc, 1) * ei_gpts';

    hem.doquad = @(hem, kernel) doquad(hem, kernel);
end

function Gpts = gauss_points(N, h, low_lim)
    % determines the Gauss points for a four point quadrature rule.

    gpts(1) = (1/2) - sqrt(15+2*sqrt(30))/(2*sqrt(35));
    gpts(2) = (1/2) - sqrt(15-2*sqrt(30))/(2*sqrt(35));
    gpts(3) = (1/2) + sqrt(15-2*sqrt(30))/(2*sqrt(35));
    gpts(4) = (1/2) + sqrt(15+2*sqrt(30))/(2*sqrt(35));

    % determines the Gauss points for all N intervals.

    for gct = 1:N
        for ell = 1:4
            Gpts((gct-1)*4 + ell,1) = (gct-1)*h + gpts(ell)*h + low_lim;
        end
    end

end

function w = gauss_weights(N, h)
    % determine the Gauss weights for a four point quadrature rule

    weights = zeros(4,1);
    weights(1) = 49*h/(12*(18 + sqrt(30)));
    weights(2) = 49*h/(12*(18 - sqrt(30)));
    weights(3) = weights(2);
    weights(4) = weights(1);

    % copy the weights to form a vector for all N intervals

    Weights = weights;
    for gct = 1:N-1
        Weights = [Weights; weights];
    end

    w = Weights;
end

function v = doquad(hem, kernel)
    % Use the quadrature rule defined in hem to evaluate the integral with
    % the kernel 'kernel' (with distribution defined in hem included).

    v = hem.distc * kernel * hem.disti;
end
