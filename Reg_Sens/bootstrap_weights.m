function weights_id = bootstrap_weights(id, varargin)
% --------------------------------------------------------------------------
% Draw multinomial or Bayesian clustered bootstrap weights
% Weights are equal within panel units
%
% Inputs:
% - id:                     obs x 1 of panel unit IDs
% - varargin (optional):    string specifying 'bayes' for Bayesian bootstrap
% Outputs:
% - weights_id:             obs x 1 bootstrap weight of each observation
%
% This version: August 27, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% --------------------------------------------------------------------------


[~,~,id_unique] = unique(id);   % Unit IDs from 1 to N
N = max(id_unique);             % Number of units

% Draw weights
if ~isempty(varargin) && strcmp(varargin{1}, 'bayes')
    weights = exprnd(1,N,1); % Bayesian bootstrap
else
    weights = mnrnd(N,ones(1,N)/N)'; % Multinomial bootstrap
end

% Weights for each observation
weights_id = weights(id_unique);

end