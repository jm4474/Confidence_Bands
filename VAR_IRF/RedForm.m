function redf = RedForm(Y_raw,tau)
% -------------------------------------------------------------------------
% Reduced-form estimators of a VAR(tau) model  
% 
% Inputs:
% - Y_raw:  T x d matrix containing the time series
% - tau:    number of lags in the VAR model
% Outputs:
% - redf:   struct with the fields listed below
% redf fields:
% - nu:     VAR intercepts
% - AL:     VAR model coefficients
% - Psi:    covariance matrix of VAR model residuals
% - eta:    VAR model residuals
% - X:      VAR model regressors
% - Y:      VAR model outcomes (after truncating tau lags)
%
% The estimation always includes a constant
%
% This version: August 24, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% -------------------------------------------------------------------------


redf = struct;


%% Definitions

redf.Y_init = Y_raw(1:tau,:); % Pre-sample outcomes
redf.Y= Y_raw((tau+1):end,:); %The rows of this matrix are Y_t'

aux = lagmatrix(Y_raw,1:1:tau);
redf.X= [ones(size(redf.Y,1),1), aux((tau+1):end,:)]; %The rows of this matrix are [1,X_{t}']


%% Least-squares regression

slopeparameters=(redf.X\redf.Y)'; %contains the d x 1 constant vector and AL


%% Generate nu and the vec(A(L)) estimators

redf.AL=slopeparameters(:,2:end);    %d x d*tau
redf.nu=slopeparameters(:,1);        %d x 1


%% Innovations covariance matrix

redf.eta = redf.Y-redf.X*slopeparameters';
redf.Psi = (redf.eta'*redf.eta)/(size(redf.eta,1));


end

