clear;
addpath('../../SimInferenceClass','../../VAR_IRF');

% Confidence bands for IRFs from Gertler & Karadi (AEJ Macro 2015)
% External IV identification

% This version: August 15, 2017
% J. L. Montiel Olea & M. Plagborg-Moller

%% Load overall settings and data

gk_settings;               % Script to prepare settings

gk_VARdata;                % Script to prepare VAR data

gk_ivdata;                 % Script to load external instrument

%% IV-specific settings

% Plug-in, Bayesian, and bootstrap inference
numdraws_supt        = 1e5; % Number of normal draws used to compute plug-in sup-t crit. val.

numdraws_boot        = 1e4; % Set to 0 if bootstrap inference undesired

verbose = true;             % Print progress of bootstrap procedure

rng(20170114);              % Seed for random number generator

% Bands

band_list            = {'Pwise',...
                        'supt',...
                        'Sidak',...
                        'Bonferroni'}; %, 'thetaproj', 'muproj'}; 
                            % Cell array of bands to be plotted, 
                            % can be any combination of: 'Pwise', 'supt', 'Sidak', 'Bonferroni', 'thetaproj', 'muproj'

legend_bands         = {'Pointwise',...
                        'Sup-t: plug-in',...
                        'Sidak',...
                        'Bonferroni'}; %, '\theta-projection', '\mu-projection'}; 
                            % Legend for bands

linestyle_supt_bands = {'-', '--'};

%% Reduced-form VAR estimation

redf                 = iv_RedForm(Y, Z, tau); 
                                      % Struct redf contains reduced-form VAR
                                      % and IV objects

%% Plug-in bands

% IV estimation of IRFs and delta method variance

[Theta, Sigmahat, p] = iv_estim(redf, hori);

% Confidence bands
sel                  = select_IRF(d, 1, hori,...
                                  plot_var, 1); 
                                      % Selection vector for IRF of interest

bands_plugin         = SimInference.bands_plugin(Theta(sel)',...
                       Sigmahat(sel,sel), p, band_list, numdraws_supt,...
                       1-signif_level);% Collection of plug-in bands

% Plot plug-in bands
figure('Unit', 'normalize', 'Position', [0.2 0.2 0.6 0.6], 'Name', 'Bands: IV, plug-in');

plot_compare({Theta(sel)'}, bands_plugin, plot_band_xlabel, plot_band_ylabel, plot_xticks, legend_bands);

%% Bootstrap band

if numdraws_boot > 0
   
 tic;
    
 vecTheta_draws_boot = iv_boot(redf, hori,...
                          numdraws_boot,...
                          verbose); toc  % Run bootstrap
    
    calibband_boot   = SimInference.calibrated_Rbands(...
                       vecTheta_draws_boot(:,sel),...
                       1-signif_level); % Calibrated bootstrap
end

%% Plot comparing sup-t bands

if numdraws_boot > 0
    
    figure('Unit', 'normalize', 'Position', [0.2 0.2 0.6 0.6], 'Name', 'Bands: IV, sup-t');
    
    plot_compare({Theta(sel)'}, {bands_plugin{find(strcmp(band_list,'supt'),1)}, calibband_boot}, plot_band_xlabel, plot_band_ylabel, plot_xticks, {'Sup-t: plug-in', 'Sup-t: bootstrap'}, linestyle_supt_bands);

end
