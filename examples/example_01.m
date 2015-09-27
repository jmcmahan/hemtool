% This example uses the default parameters for HEMtool. This results in
% a 90 degree switching forward model which computes the strain for a 
% particular PZT actuator resulting from an applied electric field input. 
% The defaults are in the file, 'default_param.m', in the HEMtool source.

clear all;


% 'hemtool' returns a structure containing functions implementing the model
% chosen by the parameters. Defaults are used for parameters not passed.
hem = hemtool(); 

% The defaults provide a delta_t parameter which defines the time between
% subsequent inputs. It is used by the algorithms which use pre-computed
% lookup tables to implement the update of the dipole fraction.
dt = hem.delta_t;

t = (0:dt:dt*200)';

freq = 1/10; % Hz
E = 5e5 + 5e5 * sin(freq * 2*pi*t + pi); % Sinusoidal input for this example


% Here the forward algorithm is computed for the given input (in this case,
% the strain, eps, is computed for the electric field, E). If E is a vector,
% the values are considered the input at each time step, so for this example,
% E(k) corresponds to t(k). The corresponding output is saved in eps. Note
% that 'hem' is also updated, as well as passed as an input. This is necessary
% if any additional calls to hem.forward are made, since the state information 
% which reflects the history of the input is stored in 'hem'.

[eps, hem] = hem.forward(hem, E);


% Here we initialize the inverse algorithm using the defaults. A separate 
% structure is needed, since the inverse algorithm has its own state variables
% to update.

% This tells the initialization routine to setup the inverse function and use
% the default material and inverse parameters. The structure implementing
% functions for the inverse will be stored in 'heminv'.

param.function_type = 'inverse';
heminv = hemtool(param);

% Now we pass the output computed above as input to the inverse algorithm to
% try to approximately recover 'E'.

[Ealt, heminv] = heminv.invert(heminv, eps);

% Finally, we'll use one more hemtool structure to compute the strain using
% the values the input field computed by the inverse, to evaluate the 
% accuracy of the inverse function for this problem.
hemalt = hemtool();
[epsalt, hemalt] = hemalt.forward(hemalt, Ealt);


% Plot input field, E, versus the input field, E2, recovered by the 
% inverse algorithm.
figure(1);
plot(t, (E-Ealt)/1e6);
legend('E-Ealt');
title('Input Error From Inverse Algorithm');
ylabel('Electric Field (MV/m)');
xlabel('Time (s)');

% Plot strain versus input field for the original input along with the strain
% computed from the input that came from the inverse algorithm.
figure(2);
plot(E/1e6, eps*100, Ealt/1e6, epsalt*100, '--');
legend('eps', 'epsalt');
title('Strain versus Electric Field');
ylabel('Strain (%)');
xlabel('Electric Field (MV/m)');


