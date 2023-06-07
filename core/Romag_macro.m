%% Romag macro

% Channels Correlation
run('channel_map.m');

% Geometry and material properties
run('constants.m');

% Heat loss estimate
run('p_loss_corr.m');

% Material Properties
[rho,alpha,Cp,k,sigma,eta,kappa,nu] = Mat_props(mean(T_mean),2);

% Scales & Dimensionless parameters
run('scale_param.m');

% Frequency estimates
run('freq_predict.m');

% Color map
run('colors.m');