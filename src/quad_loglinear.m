function hem = quad_loglinear(hem, param)
    % Calculates the grid points and quadrature weights using a loglinear
    % quadrature sampling.

    % Nc and Ni should be divisible by 2.
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
    ec_gpts_p = get_loglinear(ec_bar, ec_int, Nc/2, 1);
    ec_gpts_n = get_loglinear(0, ec_bar, Nc/2, -1);
    ec_gpts = [ec_gpts_n(1:end-1) + ec_gpts_n(2:end);
               ec_gpts_p(1:end-1) + ec_gpts_p(2:end)] / 2;

    wc = diff([ec_gpts_n(1:end-1); ec_gpts_p]);

    % Setup interaction field grid points
    ei_gpts_p = get_loglinear(0, ei_int, Ni/2, 1);
    ei_gpts_n = -flipud(ei_gpts_p);
    ei_gpts = [ei_gpts_n(1:end-1) + ei_gpts_n(2:end);
               ei_gpts_p(1:end-1) + ei_gpts_p(2:end)] / 2;

    wi = diff([ei_gpts_n(1:end-1); ei_gpts_p]);

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

function gpts = get_loglinear(lb, ub, N, K)
    gpts = zeros(N, 1);
    b = K*4.0;
    y = 0:1/N:1;
    y = y';
    a = (ub - lb) / (exp(b) - 1);
    c = (ub - exp(b)*lb) / (ub - lb);
    gpts = a*(exp(b*y) - c);
end

function v = doquad(hem, kernel)
    % Use the quadrature rule defined in hem to evaluate the integral with
    % the kernel 'kernel' (with distribution defined in hem included).

    v = hem.distc * kernel * hem.disti;
end

