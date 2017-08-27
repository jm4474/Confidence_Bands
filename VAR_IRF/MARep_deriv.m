function [Phi, dvecPhi_dvecAL] = MARep_deriv(AL, hori)
% ------------------------------------------------------
% Reduced-form IRFs and their derivative wrt. A(L) coefficients
% Uses notation of Lütkepohl (2005), chapter 3.7
%
% Inputs:
% - AL:     VAR model coefficients [A_1,...,A_p]
% - hori:   largest response horizon
% Outputs:
% - Phi:    reduced-form IRFs [Phi_1,...,Phi_hori]
% - dvecPhi_dvecAL: dvec(Phi)/dvec(AL)
%
% This version: August 24, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% ------------------------------------------------------

%% Auxiliary matrices

d = size(AL,1);
tau = size(AL,2)/d;


%% Reduced-form MA coefficients

Phi = MARep(AL, hori);


%% Compute derivative, if desired

if nargout>1
    
    % Assemble VAR companion matrix
    if tau==1
        A = AL;
    else
        A = [AL; eye(d*(tau-1)) zeros(d*(tau-1),d)];
    end

    %% Compute derivative iteratively (cf. Lütkepohl)
    deriv = zeros(d^2, d^2*tau, hori);
    for i=1:hori
       for m=0:i-1
          if m==0
              Phi_m = eye(d);
          else
              Phi_m = Phi(:,d*(m-1)+1:d*m);
          end
          A_pow = (A')^(i-1-m);
          deriv(:,:,i) = deriv(:,:,i) + kron(A_pow(1:d,:), Phi_m);
       end
    end

    %% Return vec
    dvecPhi_dvecAL = reshape(permute(deriv, [1 3 2]), d^2*hori, d^2*tau);
end

end

