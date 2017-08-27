function V = cov_clust(Z, id)
% ---------------------------------------------
% Compute cluster (sandwich) covariance matrix
%
% Inputs:
% - Z:      obs x m data matrix
% - id:     obs x 1 vector of unit IDs
% Outputs:
% - V:      m x m cluster covariance matrix
%
% This version: August 27, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% ---------------------------------------------

sums_id = zeros(max(id),size(Z,2));         % Will contain the unit sums for each column of Z
for j=1:size(Z,2) % For each column of Z...
    sums_id(:,j) = accumarray(id, Z(:,j));  % Sum by unit ID
end
V = sums_id'*sums_id;                       % Cluster sandwich matrix

end