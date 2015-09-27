% This is an example showing an INCORRECT and CORRECT way to compute the
% polarization when using HEMtool.


% On initialization, HEMtool returns a structure containing state information
% needed for computing the desired model. It is important that this state
% information be properly updated for the output to give the desired result.
% HEMtool's functions automatically update the state information and store it
% back in the model structure which is returned, so it is necessary to update
% the model structure.
%
% The following will demonstrate, using default choices for the input
% parameters.

% Input setup

t = 0:0.1:20'; % Time vector - we've set it up to have the same delta_t
               % as the default in HEMtool.
E = 5e5*sin(0:0.1:20)'; % Input field vector

% For the first example, we'll update incorrectly.

% Setup model structure, hwrong
hwrong = hemtool();


%%%%%%%%%%%%%%%%%%%%%%%%
%%% INCORRECT METHOD %%%
%%%%%%%%%%%%%%%%%%%%%%%%

epswrong = zeros(size(E));
for i = 1:length(E)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% HERE IS THE MISTAKE %%%  
    epswrong(i) = hwrong.forward(hwrong, E(i));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% EXPLANATION %%%
    % In the above, we are only saving the strain, 'eps', rather than
    % using [eps, hwrong] = ... to update the model structure, hwrong. The 
    % effect is that the dipole fraction, hwrong.xp and hwrong.xm, does not 
    % get updated, so rather than obtaining various values of 'eps' as time
    % progresses, we obtain 'eps' at this specific time for the various values
    % of E(i). This will not normally be what you want.
end


% Now we'll do it the correct way.


% Setup model structure, hright
hright = hemtool();

epsright = zeros(size(E));
for i = 1:length(E)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% HERE IS THE MISTAKE %%%  
    [epsright(i), hright] = hright.forward(hright, E(i));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% EXPLANATION %%%
    % In the above, since 'hright' is being updated, the updated values for
    % the dipole fraction, hright.xp and hright.xm, are being updated. Thus
    % the state is progressing in time and we obtain the values of "eps" for
    % various time steps.
    %
end

% Note in the above examples, where you have some known input values, 
% E(n), corresponding to times, t(n), you can simply pass the input vector
% rather than manually looping. In cases where the input depends on the
% output in some way, though, the above points are important. For example,
% if we had
% E(n) = 5*eps(n-1);
% there would be no way to compute the vector E(n) ahead of time.





% Plot for the incorrect update, fixed in time. Note the lack of hysteresis
figure(1);
plot(E/1e6, epswrong*100, 'r');
legend('incorrect');
title('Incorrect Update, Time fixed at t=t(1)');
xlabel('Electric Field (MV/m)');
ylabel('Strain (%)');

% Note for the correct update, we see hysteresis 
figure(2);
plot(E/1e6, epsright*100, 'g');
legend('correct');
title('Correct Update, Time fixed at t=t(1)');
xlabel('Electric Field (MV/m)');
ylabel('Strain (%)');

% This 3D plot will illustrate what is being computed in each case. Note
% how the curve which comes from 'epswrong' is fixed at the initial time.
% It's instructive to rotate the resulting plot to get a clear picture of
% the shape in these dimensions.

figure(3);
plot3(E/1e6, t, 100*epsright, 'g', ...
      E/1e6, t(1)*ones(size(t)), 100*epswrong, 'r');
legend('correct', 'incorrect');
title('Strain versus both Input and Time');
xlabel('Electric Field (MV/m)');
ylabel('Time (s)');
zlabel('Strain (%)');
