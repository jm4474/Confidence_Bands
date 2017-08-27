function [H, dH_dvechPsi_gamma] = iv_H_deriv(gamma, Psi)
% ------------------------------------------------------
% Compute external IV impulse vector, and its derivative
% Impulse vector is defined as
% H = (1/sqrt(gamma'*Psi^{-1}*gamma)))*gamma
%
% Inputs:
% - gamma:    d x 1
% - Psi:      d x d
% Outputs:
% - H:                      d times 1
% - dH_dvechSigma_gamma:    deriv. of H wrt. [vech(Psi); gamma]
%
% This version: August 23, 2017
% J. L. Montiel Olea and M. Plagborg-Moller
% ------------------------------------------------------


%% Calculation of H

d = length(gamma);
Psiinv_gamma = Psi\gamma;
denom = sqrt(gamma'*Psiinv_gamma);
H = gamma/denom;


%% Calculation of derivative via chain rule

% Derivative wrt. gamma
dH_ddenom = -gamma/denom^2; % "denom" is sqrt(gamma'*Psi^{-1}*gamma)
ddenom_dradicand = 1/(2*denom); % "radicand" is gamma'*Psi^{-1}*gamma
dH_dgamma = eye(d)/denom + 2*dH_ddenom*ddenom_dradicand*Psiinv_gamma';

% Derivative wrt. vec(Psi)
dradicand_dPsi = -Psiinv_gamma*Psiinv_gamma';
dH1_dvecPsi = dH_ddenom*ddenom_dradicand*dradicand_dPsi(:)'; % dH/dvec(Psi)

% Derivative wrt. [vech(Psi); gamma]
select_vech = find(tril(ones(d))); % Vector that selects linear indices of vech(Psi) from indices of vec(Psi)
mult = 2*ones(d) - diag(ones(1,d)); % Take into account that vech(Psi) influences both lower and upper triangle of Psi
dH_dvechPsi_gamma = [dH1_dvecPsi(:,select_vech).*mult(select_vech)' dH_dgamma]; % dH/d[vech(Psi); gamma]

end