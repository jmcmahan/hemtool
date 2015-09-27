% This example computes polarization for a given input field using the
% 180 degree switching model for four different cases:
%
% CASE 1: normal / lognormal distributions and thermal relaxation
% CASE 2: normal / lognormal distributions and negligible thermal relaxation
% CASE 3: general distribution and thermal relaxation
% CASE 4: general distribution and thermal negligible relaxation
%
% Note that the "general distribution" means that the distribution is 
% specified from given data. The general distributions given here do not
% necessarily correspond to the distributions approximated by normal and
% lognormal distributions, so the differences in the output should not
% be taken as representing an accurate comparison of the methods.

clear all;


% Setup the input field

E_1 = linspace(0,.5,26)';
E_2 = linspace(.5,-.5,51)';
E_3 = linspace(-.5,2,126)';
E_4 = linspace(2,-2,201)';
E_5 = linspace(-2,.4,121)';
E_6 = linspace(.4,-2,121)';
E_7 = linspace(-2,.6,131)';
E_8 = linspace(.6,-.3,46)';
E_9 = linspace(-.3,2,116)';
E_10 = linspace(2,-.7,136)';
E_11 = linspace(-.7,1,86)';
E_12 = linspace(1,-.5,76)';
E_13 = linspace(-.5,.5,51)';
E_14 = linspace(.5,-1.2,86)';
E_15 = linspace(-1.2,-.6,31)';
E_16 = linspace(-.6,-2,71)';


E = [E_1(1:end-1); E_2(1:end-1); E_3(1:end-1); E_4(1:end-1);...
     E_5(1:end-1); E_6(1:end-1); E_7(1:end-1); E_8(1:end-1);...
     E_9(1:end-1); E_10(1:end-1); E_11(1:end-1); E_12(1:end-1);...
     E_13(1:end-1); E_14(1:end-1); E_15(1:end-1); E_16(1:end)]*1e6;

clear E_1 E_2 E_3 E_4 E_5 E_6 E_7 E_8 E_9 E_10 E_11 E_12 E_13 E_14 E_15 E_16;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CASE 1:                                              %%%
%%% LOGNORMAL / NORMAL DISTRIBUTION + THERMAL RELAXATION %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This says to use the evolution equation corresponding to the 180 degree
% switching model incorporating thermal relaxation (that is, update the 
% positive dipole fraction, which will be stored in hem.xp, according to the
% rules given in the 180 degree switching model). The code in 
% 'evol_180_relax.m' gets used.
param.evolution = 'evol_180_relax_notab';
% This says to evaluate the state for a given input to give the polarization
% using the function set up in 'eval_180_polarization.m'.
param.evaluation = 'eval_180_polarization';

% Use this time step when computing lookup tables used in updating the 
% dipole fraction.
param.delta_t = 5e-2;


% These parameters derived from the code in Chapter 2 of the Smart Materials
% book

% Parameters for the normal / lognormal distributions
param.ec_bar = 7.5810023029014096e5;
param.sigma_c = 2.3946655490458121e-01*sqrt(2);
param.sigma_i = sqrt(2.9849260411205185e+10/2);

param.eta = 1.0308396819793936e+07;

% These parameters are set to 1, indicating we're using just 1 of each 
% of the normal / lognormal distributions, rather than a weighted sum.
param.sigma_c_scale = [1];  
param.mean_scale = [1];
param.sigma_i_scale = [1];
param.dist_c_weights = [1];
param.dist_i_weights = [1];

% Remanence polarization
param.P_r = 0.3;

% This value currently has no actual function, but creates a field in the
% resulting structure identifying the type of output.
param.output_type = 'polarization';

hem = hemtool(param);

% Compute the polarization
[P1, hem] = hem.forward(hem, E);


figure(1);
plot(E/1e6, P1);
legend('polarization');
title('Polarization versue Electric Field, lognormal, relaxation on');
xlabel('Electric Field (MV/m)')
ylabel('Polarization (C/m^2)')

