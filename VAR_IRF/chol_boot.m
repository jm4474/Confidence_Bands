function vecTheta_draws = chol_boot(redf, hori, numdraws, verbose)
% -------------------------------------------------------------
% Generate (homoskedastic) residual bootstrap VAR IRF estimates
% Cholesky identification
%
% Inputs:
% - redf:       struct with reduced-form VAR objects
% - hori:       largest response horizon
% - numdraws:   number of bootstrap draws
% - verbose:    true if progress should be printed on screen
% Outputs:
% - vecTheta_draws:     draws of vec(Theta) (# rows: numdraws)
%
% This version: August 24, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% -------------------------------------------------------------


%% Dimensions

d = size(redf.AL,1);
tau = size(redf.AL,2)/d;
T = size(redf.eta,1);


%% Bootstrap

vecTheta_draws = zeros(numdraws, d^2*(hori+1));
if verbose
    disp('Bootstrapping...');
end

for j=1:numdraws
      
   % Resample VAR residuals
   eta_b = datasample(redf.eta, T);
   
   % Bootstrap VAR data
   Y_b = zeros(T, d);
   Y_b(1:tau,:) = redf.Y_init;
   for t=tau+1:T
      y_lag = Y_b(t-1:-1:t-tau,:)';
      Y_b(t,:) = [1 y_lag(:)']*[redf.nu redf.AL]' + eta_b(t,:);
   end
   
   % Reduced-form calculations
   redf_b = RedForm(Y_b, tau);
   
   % Estimate IRF
   Theta_draw = chol_estim(redf_b, hori);
   vecTheta_draws(j,:) = Theta_draw(:)';
   
   % Print progress
   if verbose && mod(j, ceil(numdraws/20))==0
      fprintf('%3d%s\n',  ceil(100*j/numdraws), '%');
   end
   
end


end