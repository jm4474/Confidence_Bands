%% Data and specification settings

% File
filename = '../Data/gravity.csv'; % CSV data file

% Parameter of interest
hori = 40; % Number of years after independence

% Specifications for sensitivity analysis
controls = cell(1,5); % Cell array of names of control variables in different specifications
controls{1} = {'col_hist', 'col_a'};
controls{2} = [controls{1}, {'pop_o_log', 'pop_d_log', 'gdpcap_o_log', 'gdpcap_d_log', 'distw_log', 'contig'}];
controls{3} = [controls{2}, {'comlang_off', 'comleg'}];
controls{4} = [controls{3}, {'comcur'}];
controls{5} = [controls{4}, {'rta', 'gatt_both', 'acp_to_eu'}];


%% Load data

dat = readtable(filename); % Data table

[~,~,id] = unique(dat.dyad_id);                 % Unit IDs
Y = dat.flow_roundlog;                          % Outcome
yrsindepcap_fe = dummy_var(dat.yrsindepcap);    % dummy variables: years since independence
year_fe = dummy_var(dat.year);                  % dummy variables: year
X = yrsindepcap_fe(:,hori);                     % Covariate of interest

obs = length(id); % Total number of observations


%% Control variable data matrix

k = length(controls); % Number of specifications

% Generate list of all (unique) control variables
control_vars_all = cell(0);
for i=1:k % For each specification...
    control_vars_all = [control_vars_all, controls{i}];
end
control_vars_unique = unique(control_vars_all);         % List of all unique control variable names
num_control_vars_unique = length(control_vars_unique);  % Number of unique control variables

% Data matrix with all control variables
W_all = [dat{:,control_vars_unique} ones(obs,1) yrsindepcap_fe(:,[1:hori-1 hori+1:size(yrsindepcap_fe,2)]) year_fe]; % Covariates, constant, fixed effects
num_controls = size(W_all,2);                               % Number of regressors
num_dummy_controls = num_controls-num_control_vars_unique;  % Number of dummy variable regressors


%% List of specifications for sensitivity analysis

% Create logical matrix indicating which controls are used in each specification
controls_ind = false(k, num_controls);
for i=1:k % For each specification...
    controls_ind(i,num_control_vars_unique+1:end) = true;                   % Always include dummy variables
    for l=1:length(controls{i}) % For every control variable in specification...
        controls_ind(i,strcmp(control_vars_unique, controls{i}{l})) = true; % Also include relevant control
    end
end
