function [Theta, Sigmahat, p] = chol_estim(redf, hori)
% ------------------------------------------------------------------
% Estimate VAR IRFs by Cholesky and compute delta method variance
%
% Inputs:
% - Y:      T x d VAR data
% - Z:      T x 1 IV data (some entries may be missing)
% - tau:    number of VAR lags
% - hori:   largest IRF horizon to compute
% Outputs:
% - Theta:      structural IRFs [Theta_0, ..., Theta_hori]
% - Sigmahat:   delta method var-cov matrix of vec(Theta)
%
% This version: August 23, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% ------------------------------------------------------------------


%% Impact impulse responses

[H, dH_dvechPsi] = chol_H_deriv(redf.Psi); % Impact responses H and derivative


if nargout == 1
    
    %% Structural IRFs
    
    % Compute only structural IRFs, if derivatives not desired
    Theta = IRF_deriv(redf.AL, hori, H, [], []);

else
    
    %% Asymptotic variance via delta method
    
    % In the following we use the notation
    % mu = [vec(AL); vech(Psi)];
    
    % Dimensions of components of mu
    d = size(redf.AL,1);
    tau = size(redf.AL,2)/d;
    dim_vecAL = d^2*tau;
    dim_vechPsi = d*(d+1)/2;
    
    % Derivative of Phi and H wrt. mu
    [~, dvecPhi_dvecAL] = MARep_deriv(redf.AL, hori); % dvec(Phi)/dvec(AL)
    dvecPhi_dmu = [dvecPhi_dvecAL, zeros(size(dvecPhi_dvecAL,1),dim_vechPsi)]; % dvec(Phi)/dmu
    dvecH_dmu = [zeros(size(dH_dvechPsi,1),dim_vecAL), dH_dvechPsi]; %dvec(H)/dmu
    
    % Derivative of structural IRFs wrt. mu
    [Theta, dvecTheta_dmu] = IRF_deriv(redf.AL, hori, H, dvecPhi_dmu, dvecH_dmu); % Theta and dvec(Theta)/dmu
    
    % Variance of mu
    Omegahat = chol_Var_mu(redf.X,redf.eta);
    p = size(Omegahat,1); % Number of reduced-form parameters
    
    % Delta method variance of structural IRFs
    Sigmahat = dvecTheta_dmu*Omegahat*dvecTheta_dmu'; % Variance of IRFs

end

end