function V = variance_vec(Y, X)
% ----------------------------------------------
% Returns sample var-cov matrix of vec(Y_i*X_i')
%
% Inputs:
% Y: n x p matrix [Y_1, ..., Y_n]'
% X: n x k matrix [X_1, ..., X_n]'
% Outputs:
% V: pk x pk covariance matrix of vec(Y_i*X_i')
% ----------------------------------------------


%% Compute covariance matrix
R = repmat(Y, 1, size(X,2)) .* kron(X, ones(1, size(Y,2)));
V = cov(R, 1);

end