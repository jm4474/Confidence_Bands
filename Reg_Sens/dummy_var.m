function [dummies,N] = dummy_var(vari)
% ---------------------------------------------------------
% Generate dummy variables for categorical variable
%
% Inputs:
% - vari:      obs x 1 data vector for categorical variable
% Outputs:
% - dummies:   obs x (N-1) matrix of dummies (dropping one)
% - N:         number of unique values in vari
%
% This version: August 27, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% ---------------------------------------------------------

unique_val = unique(vari);  % Unique values
N = length(unique_val);     % Number of unique values

% Create one dummy column for every unique value, but drop one
dummies = zeros(length(vari), N-1);
for i=2:N
    dummies(:,i-1) = (vari==unique_val(i));
end

end