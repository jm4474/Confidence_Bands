function redf = iv_RedForm(Y_raw, Z_raw, tau)
% ------------------------------------------------------
% Reduced-form estimators for VAR and IV
%
% Inputs:
% - Y_raw:  T x d matrix containing the time series
% - Z_raw:  T x 1 vector containing the IV (some elements can be missing)
% - tau:    number of lags in the VAR model
% Outputs:
% - redf:   struct with the fields listed in "RedForm.m" file as well as...
% redf additional fields:
% - iv_sample:  logical vector indicating non-missing IV observations
% - beta:       coefficients for projection of IV on lagged Ys
% - v:          residuals in above projection
%
% This version: August 24, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% ------------------------------------------------------

%% Reduced-form VAR

redf = RedForm(Y_raw, tau);


%% Reduced-form IV calculations

% IV sample
Z_trunc = Z_raw(tau+1:end); % Remove pre-sample observations from IV
redf.iv_sample = ~isnan(Z_trunc); % Sample where IV is non-missing
T_iv = sum(redf.iv_sample);  % Length of IV sample

% Compute gamma = covariance of IV with reduced-form residuals
redf.gamma = redf.eta(redf.iv_sample,:)'*Z_trunc(redf.iv_sample)/T_iv;

% Compute v = residual of IV after projecting out lagged Ys
redf.beta = redf.X(redf.iv_sample,:)\Z_trunc(redf.iv_sample); % Projection of Z on X
v_sample = Z_trunc(redf.iv_sample) - redf.X(redf.iv_sample,:)*redf.beta;
redf.v = nan(size(redf.eta,1),1);
redf.v(redf.iv_sample) = v_sample;

end