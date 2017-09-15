%% Load VAR data

var_file = '../Data/VAR_data.csv';      % Data settings

var_cols = [3 4 6 11];                  % Column numbers of VAR variables
 
var_dat = importdata(var_file);         % Read from file

Y       = var_dat.data(:,var_cols);     % Define variables

d       = size(Y,2);
