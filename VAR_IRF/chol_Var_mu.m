function Omegahat = chol_Var_mu(X,eta)
% -------------------------------------------------------------------------
% Computes the variance of [vec(Ahat); vech(Psihat)]
% assuming residuals are i.i.d. homoskedastic normal
% 
% Inputs:
% - X:      T x (d*tau + 1) (d is the dimension of the VAR, tau the lag length)
% - eta:    T x d reduced-form residuals
% Outputs:
% - Omegahat: asymptotic variance of [vec(Ahat)',vech(Psihat)']'
%
% This version: August 24, 2017
% Last revised by M. Plagborg-Moller
% -------------------------------------------------------------------------


%% Variance of ALhat and vec(Psihat) on VAR sample

d = size(eta,2);
[SUR_var, vecPsihat_var] = SUR_Var(X, eta); % Variance formula for SUR regression
ALhat_var = SUR_var(d+1:end,d+1:end); % Variance of ALhat (remove intercepts)


%% Stitch together variance matrix for muhat

select_vech = vec_to_vech(d); % Vector that selects linear indices of vech(Psi) from indices of vec(Psi)

% Variance matrix of muhat
Omegahat = blkdiag(ALhat_var, ...
                   vecPsihat_var(select_vech,select_vech));


end
 
   
