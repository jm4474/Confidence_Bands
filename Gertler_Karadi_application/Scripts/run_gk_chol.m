clear;
addpath('../../SimInferenceClass','../../VAR_IRF');

% Confidence bands for IRFs from Gertler & Karadi (AEJ Macro 2015)
% Recursive identification

% This version: August 24, 2017
% J. L. Montiel Olea & M. Plagborg-Moller


%% Load overall settings and data

gk_settings;              % Script to prepare settings

gk_VARdata;               % Script to prepare VAR data


%% Cholesky-identification specific settings

% Plotting

plot_shock     = 3; % Plot this shock

% Plug-in, Bayesian, and bootstrap inference

numdraws_supt  = 1e5;    % Number of normal draws used to compute plug-in sup-t crit. val.

numdraws_boot  = 1e4;    % Set to 0 if bootstrap inference undesired

numdraws_Bayes = 1e4;    % Set to 0 if Bayesian inference undesired

verbose        = true;   % Print progress of bootstrap/Bayes procedure

rng(20170114);           % Seed for random number generator

% Bands

band_list      = {'Pwise',...
                  'supt',...
                  'Sidak',...
                  'Bonferroni'}; %, 'thetaproj', 'muproj'};
                                 % Cell array of bands to be plotted, 
                                 % can be any combination of: 'Pwise', 'supt', 'Sidak', 'Bonferroni', 'thetaproj', 'muproj'
                                 
legend_bands   = {'Pointwise',...
                  'Sup-t: plug-in',...
                  'Sidak',...
                  'Bonferroni'}; %, '\theta-projection', '\mu-projection'}; 
                                 % Legend for bands
                                 
linestyle_supt_bands = {'-', '--', '-.'};

%% Reduced-form VAR estimation

redf                 = RedForm(Y, tau); % Struct redf contains reduced-form VAR objects
                                        % Y  :  vector of VAR variables
                                        % tau:  number of VAR lags
%% Plug-in bands

% Cholesky estimation of IRFs and delta method variance

[Theta, Sigmahat, p] = chol_estim(redf, hori);

% Confidence bands

sel                  = select_IRF(d, d, hori, ...
                                 plot_var, plot_shock); 
                                        % Selection vector for IRF of interest
                                        
bands_plugin         = SimInference.bands_plugin(Theta(sel)',...
                       Sigmahat(sel,sel), p, band_list, numdraws_supt,...
                       1-signif_level); % Collection of plug-in bands

% Plot plug-in bands

figure('Unit', 'normalize', 'Position', [0.2 0.2 0.6 0.6], 'Name', 'Bands: Cholesky, plug-in');

plot_compare({Theta(sel)'}, bands_plugin, plot_band_xlabel, plot_band_ylabel, plot_xticks, legend_bands);

% Store sup-t band for later plotting

bands_supt_band      = {bands_plugin{find(strcmp(band_list,'supt'),1)}};

bands_supt_legend    = {'Sup-t: plug-in'};

%% Bootstrap band

if numdraws_boot > 0

 tic;
 
 vecTheta_draws_boot = chol_boot(redf, hori,...
                          numdraws_boot,...
                          verbose); toc            % Run bootstrap
                                                   % Standard Homoskedastic 
                                                   % residual bootstrap
 calibband_boot      = SimInference.calibrated_Rbands(...
                          vecTheta_draws_boot(:,sel),...
                          1-signif_level);         % Calibrated bootstrap

 %Store band for later plotting
    
 bands_supt_band     = [bands_supt_band {calibband_boot}];
    
 bands_supt_legend   = [bands_supt_legend {'Sup-t: bootstrap'}];
    
end

%% Bayes band

if numdraws_Bayes > 0

tic;
vecTheta_draws_Bayes = chol_Bayes(redf, hori, numdraws_Bayes, verbose); toc;
                                                   % Run bootstrap
                                                   % Uhlig's posterior
                                                   % draws

calibband_Bayes      = SimInference.calibrated_Rbands(...
                          vecTheta_draws_Bayes(:,sel),...
                          1-signif_level);         % Calibrated Bayes

 % Store band for later plotting
 
 bands_supt_band     = [bands_supt_band {calibband_Bayes}];
 
 bands_supt_legend   = [bands_supt_legend {'Sup-t: Bayes'}];
    
end

%% Plot comparing sup-t bands

if numdraws_boot + numdraws_Bayes > 0

    figure('Unit', 'normalize', 'Position', [0.2 0.2 0.6 0.6], 'Name', 'Bands: IV, sup-t');
    
    plot_compare({Theta(sel)'}, bands_supt_band, plot_band_xlabel, plot_band_ylabel, plot_xticks, bands_supt_legend, linestyle_supt_bands(1:length(bands_supt_band)));

end
