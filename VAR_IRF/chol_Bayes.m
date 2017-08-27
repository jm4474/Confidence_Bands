function vecTheta_draws = chol_Bayes(redf, hori, numdraws, verbose)
% ---------------------------------------------------------------
% Bayesian inference on VAR IRFs by Cholesky
%
% Inputs:
% - redf:       struct with reduced-form VAR objects
% - hori:       largest response horizon
% - numdraws:   number of Bayes draws
% - verbose:    true if progress should be printed on screen
% Outputs:
% - vecTheta_draws:     draws of vec(Theta) (# rows: numdraws)
%
% This version: August 24, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% ---------------------------------------------------------------


%% Reduced-form posterior draws

nuAL = [redf.nu, redf.AL];
PosteriorDraws = UhligPosterior(nuAL(:),redf.Psi(:),redf.X,numdraws);


%% Structural IRF draws

d = size(redf.AL,1);
vecTheta_draws = zeros(numdraws, d^2*(hori+1));
if verbose
    disp('Transforming Bayes draws...');
end

for j=1:numdraws % Cycle through reduced-form draws
   
   % Draw of H matrix
   Psi_draw = PosteriorDraws.Psi(:,:,j);
   H_draw = chol_H_deriv(Psi_draw);
   
   % Draw of structural IRFs
   Theta_draw = IRF_deriv(PosteriorDraws.A(:,:,j), hori, H_draw, [], []);
   vecTheta_draws(j,:) = Theta_draw(:)';
   
   % Print progress
   if verbose && mod(j, ceil(numdraws/20))==0
      fprintf('%3d%s\n',  ceil(100*j/numdraws), '%');
   end
   
end


end