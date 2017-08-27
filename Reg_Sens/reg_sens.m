function [beta_hats, Sigma_hat] = reg_sens(Y, X, W, id, controls_ind)
% ---------------------------------------------
% Regression sensitivity analysis
%
% Inputs:
% - Y:              obs x 1 vector of outcomes
% - X:              obs x 1 vector of covariate of interest
% - W:              obs x m data matrix of all controls
% - id:             obs x 1 list of panel unit IDs
% - controls_ind:   k x m logical matrix indicating which controls are used in
%                   the k different specifications
% Outputs:
% - beta_hats:      k x 1 vector of coefficient of interest across specifications
% - Sigma_hat:      k x k delta method var-cov matrix of beta_hats
%
% This version: August 27, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% ---------------------------------------------


k = size(controls_ind,1);

beta_hats = zeros(k,1);
x_tilde_sqs = zeros(k,1);
scores = zeros(length(Y),k);


%% Execute each regression specification

for j=1:k % For ever specification...
    
    if nargout>1 % If var-cov matrix is desired
        
        % Residual regression, using the list of controls for this specification
        [the_beta_hat, the_X_tilde, the_resid] = resid_reg(Y, X, W(:,controls_ind(j,:)));
        
        x_tilde_sqs(j) = the_X_tilde'*the_X_tilde;  % Store X'X (residualized)
        scores(:,j) = the_X_tilde.*the_resid;       % Store regression scores
        
    else % If only point estimate is required
        
        % Regression
        the_beta_hat = resid_reg(Y, X, W(:,controls_ind(j,:)));
        
    end
    
    beta_hats(j) = the_beta_hat; % Point estimate of interest
    
end


%% Compute grand var-cov matrix Sigma_hat, if desired

if nargout>1
    
    % Cluster sandwich matrix
    sandwich = cov_clust(scores, id);
    
    % Var-cov matrix of beta_hats
    Sigma_hat = bsxfun(@rdivide, bsxfun(@ldivide, x_tilde_sqs, sandwich), x_tilde_sqs');
    
end

end