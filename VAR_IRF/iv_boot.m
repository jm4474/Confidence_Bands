function vecTheta_draws = iv_boot(redf, hori, numdraws, verbose)
% -------------------------------------------------------------
% Generate (homoskedastic) residual bootstrap VAR IRF estimates
% External IV identification
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
T_iv = sum(redf.iv_sample);


%% Bootstrap

vecTheta_draws = zeros(numdraws, d*(hori+1));
if verbose
    disp('Bootstrapping...');
end

for j=1:numdraws
   
   % Initialize residuals
   eta_b = nan(T,d);
   v_b = nan(T,1);
   
   % Resample VAR and IV residuals on IV sample
   res_b = datasample([redf.eta(redf.iv_sample,:) redf.v(redf.iv_sample)], T_iv);
   eta_b(redf.iv_sample,:) = res_b(:,1:d);
   v_b(redf.iv_sample) = res_b(:,end);

   % Resample VAR residuals on remaining sample
   eta_b_remaining = datasample(redf.eta(~redf.iv_sample,:), sum(~redf.iv_sample)); % Resample VAR residuals on this sample
   eta_b(~redf.iv_sample,:) = eta_b_remaining;
   
   % Bootstrap VAR and IV data
   Y_b = zeros(T, d);
   Z_b = zeros(T, 1);
   Y_b(1:tau,:) = redf.Y_init;
   for t=tau+1:T
      y_lag = Y_b(t-1:-1:t-tau,:)';
      Y_b(t,:) = [1 y_lag(:)']*[redf.nu redf.AL]' + eta_b(t,:);
      Z_b(t) = [1 y_lag(:)']*redf.beta + v_b(t);
   end
   
   % Reduced-form calculations
   redf_b = iv_RedForm(Y_b, Z_b, tau);
   
   % Estimate IRF
   Theta_draw = iv_estim(redf_b, hori);
   vecTheta_draws(j,:) = Theta_draw(:)';
   
   % Print progress
   if verbose && mod(j, ceil(numdraws/20))==0
      fprintf('%3d%s\n',  ceil(100*j/numdraws), '%');
   end
   
end


end