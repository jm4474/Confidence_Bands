function [Theta, dvecTheta_dmu] = IRF_deriv(AL, hori, H, dvecPhi_dmu, dvecH_dmu)
% ------------------------------------------------------ 
% Structural IRFs and their derivatives
% Theta_i = Phi_i*H, i=0,1,2,... (H need not be square matrix),
% where Phi_i is matrix of reduced-form IRFs at horizon i.
% Uses notation of Lütkepohl (2005), chapter 3.7.
%
% Inputs:
% - AL: VAR model coefficients [A_1,...,A_p]
% - hori: forecast horizon
% - H: impact response matrix (may be d x m)
% - dvecPhi_dmu: dvec(Phi)/dmu, may be inputted as empty set if derivative not wanted
% - dvecH_dmu: dvec(H)/dmu, may be inputted as empty set if derivative not wanted
% Outputs:
% - Theta: structural IRFs
% - dvecTheta_dmu: dvec(Theta)/dmu
%
% This version: August 23, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% ------------------------------------------------------


%% Dimensions

d = size(AL,1); % Dimension of VAR
m = size(H,2); % Number of columns of H
p = size(dvecPhi_dmu,2); % Number of reduced-form parameters mu


%% Reduced-form IRFs

Phi = MARep(AL, hori); % Reduced-form IRFs [Phi_1,...,Phi_hori]


%% Loop over horizons to compute structural IRFs and derivatives

Theta = zeros(d, m*(hori+1)); % Will contain [Theta_0,...,Theta_hori]
dvecTheta_dmu_reshape = zeros(d*m, p, hori+1); % Will contain dvec(Theta_l)/dmu for each l=0,...,hori

% Reshape dvec(Phi)/dmu array if derivatives desired
if nargout > 1
    dvecPhi_dmu_reshape = reshape(dvecPhi_dmu, d^2, hori, p); % Reshape dvec(Phi)/dmu array
    dvecPhi_dmu_reshape = permute(dvecPhi_dmu_reshape, [1 3 2]); % Permute indices of dvec(Phi)/dmu array
end

for l=0:hori
    
    % Compute structural IRF at horizon l
    if l==0 % Impact
       Phi_l = eye(d); % Reduced-form impact responses
    else % Other horizons
       Phi_l = Phi(:,d*(l-1)+1:d*l);
    end
    Theta(:,m*l+1:m*(l+1)) = Phi_l*H; % Structural responses Theta_l
    
    % Compute derivative of structural IRF at horizon l, if desired
    if nargout > 1
        if l == 0
            dvecPhi_l_dmu = zeros(d^2,p); % Reduced-form impact responses don't depend on AL
        else
            dvecPhi_l_dmu = dvecPhi_dmu_reshape(:,:,l);
        end
        % Use chain rule to compute derivative
        dvecTheta_dmu_reshape(:,:,l+1) = kron(H',eye(d))*dvecPhi_l_dmu + kron(eye(m),Phi_l)*dvecH_dmu; % dvec(Theta_l)/dmu
    end
    
end


%% Reshape derivatives
dvecTheta_dmu = reshape(permute(dvecTheta_dmu_reshape, [1 3 2]), d*m*(hori+1), []);

end