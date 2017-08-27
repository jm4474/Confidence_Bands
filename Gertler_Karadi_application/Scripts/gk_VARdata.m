%% Load VAR data

% Data settings
var_file = '../Data/VAR_data.csv';
var_cols = [3 4 6 11]; % Column numbers of VAR variables

% Read from file
var_dat = importdata(var_file);

% Define variable
Y = var_dat.data(:,var_cols);
d = size(Y,2);
