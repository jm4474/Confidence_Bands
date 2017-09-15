clear;
addpath('../../Reg_Sens', '../../SimInferenceClass');

% Confidence bands for sensitivity analysis for Head, Mayer & Ries (JIE 2010)
% This version: August 27, 2017
% J. L. Montiel Olea & M. Plagborg-Moller


%% Settings

% Plug-in, Bayesian, and bootstrap inference
numdraws_supt = 1e5;    % Number of normal draws used to compute plug-in sup-t crit. val.
numdraws_boot = 1e2; %1e4    % Set to 0 if bootstrap inference undesired
numdraws_Bayes = 1e2; %1e4   % Set to 0 if Bayesian inference undesired
verbose = true;         % Print progress of Bayes/bootstrap procedure
rng(20170827);          % Seed for random number generator

% Plotting
signif_level = 0.1; % Significance level
control_labels = {'colonial', '+gravity', '+language/legal', '+currency', '+trade agreem'}; % Labels of specifications
bands_linestyle = {'-', '-', '--', '-.'}; % Linestyles of bands
bands_linewidth = [2 1 1 1]; % Line widths of bands


%% Load data

hmr_data; % Load data and specification settings


%% Point estimates and plug-in bands

% Point estimates and delta method var-cov matrix
[beta_hats, Sigma_hat] = reg_sens(Y, X, W_all, id, controls_ind);

% Compute plug-in bands
bands_plugin = SimInference.bands_plugin(beta_hats', Sigma_hat, [], {'Pwise', 'supt'}, numdraws_supt, 1-signif_level);

% Store bands for later plotting
bands_plot = bands_plugin;
bands_legend = {'Pointwise', 'Sup-t: plug-in'};


%% Bootstrap inference
    
if numdraws_boot>0

tic;
% Bootstrap draws
disp('Bootstrapping...');
beta_hat_draws_boot = bootstrap_reg_sens(Y, X, W_all, id, controls_ind, numdraws_boot, verbose);
toc;

% Calibrated bootstrap band
calibband_boot = SimInference.calibrated_Rbands(beta_hat_draws_boot, 1-signif_level);

% Store bands for later plotting
bands_plot = [bands_plot {calibband_boot}];
bands_legend = [bands_legend {'Sup-t: bootstrap'}];

end


%% Bayesian inference

if numdraws_Bayes>0

tic;
% Bayesian bootstrap draws
disp('Taking Bayes draws...');
beta_hat_draws_Bayes = bootstrap_reg_sens(Y, X, W_all, id, controls_ind, numdraws_Bayes, verbose, 'bayes');
toc;

% Calibrated Bayes band
calibband_Bayes = SimInference.calibrated_Rbands(beta_hat_draws_Bayes, 1-signif_level);

% Store bands for later plotting
bands_plot = [bands_plot {calibband_Bayes}];
bands_legend = [bands_legend {'Sup-t: Bayes'}];
    
end


%% Plot bands

figure('Unit', 'normalize', 'Position', [0.1 0.1 0.8 0.8], 'Name', 'Sup-t bands');

plot_compare(beta_hats, bands_plot, control_labels, bands_legend, bands_linestyle(1:length(bands_plot)), bands_linewidth(1:length(bands_plot)));

