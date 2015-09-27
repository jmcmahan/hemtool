function hem = hemtool(param)
%HEMTOOL    Setup code for implementations of the homogenized energy model.
%

    % Set up defaults for parameters not provided
    hem = struct();
    if nargin == 0
        param = '';
    end
    param = setup_parameters(param);

    quadrature = param.quadrature;
    distribution = param.distribution;
    evolution = param.evolution;
    evaluation = param.evaluation;


    % Set up the quadrature / grid points
    hem = feval(quadrature, hem, param);

    % Set up the distributions
    hem = feval(distribution, hem, param);

    % Set up the dipole evolution
    hem = feval(evolution, hem, param);

    % Set up the quadrature evaluation
    hem = feval(evaluation, hem, param);

    % Save desired parameters
    hem = save_params(hem, param);

    % Set up inverse, if applicable
    if strcmp(hem.function_type, 'inverse')
        inverse = hem.inverse;
        hem = feval(inverse, hem, param);
    end

end


function param = setup_parameters(param)
    % This function fills in any unlisted parameters with default values.
    % It also checks whether the code is being run from MATLAB or Octave
    % to address any compatibility issues.
    dflt = default_param();


    % Check inputs for problems
    if ~isstruct(param)
        disp('Warning: invalid param structure passed, using default');
        param = struct();
    end


    % Copy over fields from the default parameters if they don't exist
    % in 'param'
    fnames = fieldnames(dflt);
    for i = 1:length(fnames)
        if ~isfield(param, char(fnames(i)))
            val = getfield(dflt, char(fnames(i)));
            param = setfield(param, char(fnames(i)), val);
        end
    end

    % When making this code more full-featured, it will be good to put
    % a check here to ensure that acceptable parameter values were entered.

    param.octave = exist('OCTAVE_VERSION');
end
    

function hem = save_params(hem, param)
    hem.quadrature = param.quadrature;
    hem.distribution = param.distribution;
    hem.evolution = param.evolution;
    hem.evaluation = param.evaluation;

    hem.delta_t = param.delta_t;

    hem.output_type = param.output_type;
    hem.function_type = param.function_type;
    if strcmp(hem.function_type, 'inverse')
        hem.inverse = param.inverse;
    end

end
