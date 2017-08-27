%% Overall VAR and plotting settings

% VAR specification
tau = 12; % Number of lags in VAR
hori = 36; % Number of impulse response horizons to compute/plot

% Plot settings
signif_level = 0.32; % Pointwise/simultaneous signif. level
var_names = {'log CPI', 'log IP', '1-yr gov bond yield', 'EBP'}; % Variable names in plot
plot_var = 2; % Plot bands for this variable (and the shock chosen above)
plot_xticks = 0:12:hori; % Horizontal axis tick marks
plot_band_xlabel = 'months after shock'; % Horizontal axis label (not shown on pointwise plots)
plot_band_ylabel = 'percent'; % Vertical axis label (not shown on pointwise plots)