function Omegahat = iv_Var_mu(X,eta,v)
% -------------------------------------------------------------------------
% Computes the variance of [vec(Ahat); vech(Psihat); gammahat]
% allowing for unequal samples for Y and Z
% assuming residuals are i.i.d. homoskedastic normal
% 
% Inputs:
% - X:      T x (d*tau + 1) (d is the dimension of the VAR, tau the lag length)
% - eta:    T x d reduced-form residuals
% - v:      T x 1 
% Outputs:
% - Omegahat: asymptotic variance of [vec(Ahat)',vech(Psihat)',vec(gammahat)']'
%
% This version: August 23, 2017
% Last revised by M. Plagborg-Moller
% -------------------------------------------------------------------------


%% Variance of ALhat and vec(Psihat) on VAR sample

[T,d] = size(eta);
[SUR_var_VARsample, vecPsihat_var_VARsample] = SUR_Var(X, eta); % Variance formula for SUR regression
ALhat_var_VARsample = SUR_var_VARsample(d+1:end,d+1:end); % Variance of ALhat (remove intercepts)


%% Variance-covariance of vec(Psihat) and gammahat on IV sample

% IV sample
iv_sample = ~isnan(v);
T_iv = sum(iv_sample);

% Var-cov matrix of [vec(Psihat); gammahat]
varcov_IVsample = variance_vec(eta(iv_sample,:), [eta(iv_sample,:) v(iv_sample)])/T_iv;


%% Stitch together variance matrix for muhat

select_vech = vec_to_vech(d); % Vector that selects linear indices of vech(Psi) from indices of vec(Psi)
T_ratio = T_iv/T; % Ratio of sample sizes

% Variance matrix of muhat
mu_var = blkdiag(ALhat_var_VARsample, ...
                 [vecPsihat_var_VARsample(select_vech,select_vech), T_ratio*varcov_IVsample(select_vech,end-d+1:end);
                  T_ratio*varcov_IVsample(end-d+1:end,select_vech), varcov_IVsample(end-d+1:end,end-d+1:end)]);

% Zero out negative eigenvalues to ensure positive semidefiniteness
[V,D] = eig(mu_var);
Omegahat = V*max(D,0)*V';

end
 