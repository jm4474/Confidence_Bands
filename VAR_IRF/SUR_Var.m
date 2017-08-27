function [betahat_var, vecPsihat_var] = SUR_Var(X, res)
% ---------------------------------------------------
% Returns variance of betahat and Psihat in SUR model
% Y_i = beta*X_i + u_i, u_i ~ N(0, Psi)
% under homoskedasticity
%
% Inputs:
% - X:              T x k regressor matrix
% - res:            T x d residual matrix
% Outputs:
% - betahat_var:    variance of betahat
% - vecPsihat_var:  variance of vec(Psihat)
%
% This version: August 24, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% ---------------------------------------------------

S_X = X'*X; % Second moment of regressors
Psihat = cov(res); % Innovation variance matrix

betahat_var = kron(inv(S_X), Psihat); % Regression coefficient variance
vecPsihat_var = variance_vec(res, res)/size(res,1); % Variance of vec(Psihat)

end
