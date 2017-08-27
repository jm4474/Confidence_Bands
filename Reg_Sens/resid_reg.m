function [beta_hat, X_tilde, resid] = resid_reg(Y, X, W)
% ---------------------------------------------
% Residual regression
%
% Inputs:
% - Y:      obs x 1 outcome vector
% - X:      obs x 1 covariate of interest
% - W:      obs x k data matrix of other regressors
% Outputs:
% - beta_hat:  1 x 1 coefficient of interest
% - X_tilde:   obs x 1 residualized covariate of interest
% - resid:     obs x 1 regression residuals
%
% This version: August 27, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% ---------------------------------------------
    

%% Regression

% Remove columns of W that are all zeros
zero_cols_W = all(W==0,1);
W = W(:,~zero_cols_W);

beta_hat_full = [X W]\Y;        % Least squares
beta_hat = beta_hat_full(1);    % Coefficient of interest


%% Residualized covariate and regression residuals

if nargout > 1
    
    % Regression residuals
    resid = Y - [X W]*beta_hat_full;
    
    % Residualized covariate of interest
    gamma_hat = W\X;
    X_tilde = X - W*gamma_hat;
    
end

end