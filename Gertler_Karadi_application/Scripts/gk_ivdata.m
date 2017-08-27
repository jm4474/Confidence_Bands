%% Load IV data

% Data settings
iv_file = '../Data/factor_data.csv';
iv_col = 4; % Column number for external IV

% Read from file
iv_dat = importdata(iv_file);

% Define variable
Z = iv_dat.data(:,iv_col);
