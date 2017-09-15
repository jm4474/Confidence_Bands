%% Load IV data

iv_file = '../Data/factor_data.csv';  % Data settings

iv_col  = 4;                          % Column number for external IV

iv_dat  = importdata(iv_file);        % Read from file

Z       = iv_dat.data(:,iv_col);      % Define variable