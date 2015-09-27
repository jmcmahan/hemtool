function hem = dist_general(hem, param)
    % Construct a quadrature rule incorporating a general distribution
    % specified as a product of an interaction field distribution with
    % a coercive field distribution.

    % NOTE: Make sure the general distribution matches with the 
    % input field grid points (hem.ei_gpts and hem.ec_gpts) which are
    % generated by the quadrature rule specified.


    % These 4 parameters are specific to the general distribution
    % and required.
    ec_int = param.ec_int;
    ei_int = param.ei_int;
    % Coercive field distribution, given as a column vector
    nu_c = param.nu_c;
    % Interaction field distribution over positive, given as a 
    nu_i_pos = param.nu_i_pos;

    wi = hem.wi;
    wc = hem.wc;

    Ni = hem.Ni;
    Nc = hem.Nc;

    quadi = [];
    quadc = [];


    % Build normal interaction distribution combined with quadrature weights

    nu_i_neg = flipud(nu_i_pos);
    nu_i = [nu_i_neg; nu_i_pos];
    quadi = [quadi, nu_i .* wi];

    % Build lognormal coercive distribution combined with quadrature weights

    quadc = [quadc, nu_c .* wc];


    hem.disti = quadi;
    hem.distc = quadc';

    % Normalize the distribution to integrate to 1

    intconst = sqrt(hem.distc * ones(Nc,Ni) * hem.disti);
    hem.disti = hem.disti / intconst;
    hem.distc = hem.distc / intconst;

end
