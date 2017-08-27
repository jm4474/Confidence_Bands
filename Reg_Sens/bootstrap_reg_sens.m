function beta_hat_draws = bootstrap_reg_sens(Y, X, W_all, id, controls_ind, numdraws, verbose, varargin)
% --------------------------------------------------------------------------
% Multinomial or Bayesian bootstrap for regression sensitivity analysis
%
% Inputs:
% - Y:              obs x 1 of outcomes
% - X:              obs x 1 vector of covariate of interest
% - W_all:          obs x m matrix of all other regressors
% - id:             obs x 1 vector of panel unit IDs
% - controls_ind:   k x m logical matrix indicating which controls are used in
%                   the k different specifications
% - numdraws:       number of bootstrap draws
% - verbose:        logical, true if progress should be printed to screen
% - varargin (optional):  string specifying 'bayes' for Bayesian bootstrap
% Outputs:
% - beta_hat_draws: numdraws x k matrix of bootstrap draws
%
% This version: August 27, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% --------------------------------------------------------------------------


% Matrix of bootstrapped parameters of interest, draws along rows
beta_hat_draws = zeros(numdraws, size(controls_ind,1));

for j=1:numdraws % For each bootstrap repetition...

    weights = bootstrap_weights(id, varargin{:});   % Draw bootstrap weights
    weights_sqrt = sqrt(weights(weights>0));        % Square root of nonzero weights

    % Run sensitivity analysis on weighted data
    beta_hat_draw = reg_sens(Y(weights>0).*weights_sqrt, ...
                             X(weights>0,:).*weights_sqrt, ...
                             W_all(weights>0,:).*weights_sqrt, ...
                             id, ...
                             controls_ind);

    % Store coefficients
    beta_hat_draws(j,:) = beta_hat_draw';

    % Print progress
    if verbose && mod(j, ceil(numdraws/20))==0
      fprintf('%3d%s\n',  ceil(100*j/numdraws), '%');
    end
        
end

end