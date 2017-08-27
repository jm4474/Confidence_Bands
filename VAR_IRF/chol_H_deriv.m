function [H, dvecH_dvechPsi] = chol_H_deriv(Psi)
% ------------------------------------------------------
% Derivative of Cholesky decomposition
% Uses notation of Lütkepohl (2005), chapter 3.7
%
% Inputs:
% - Psi: reduced-form var-cov matrix
% Outputs:
% - H: chol(Psi) lower triangular, H*H'=Psi
% - dvecH_dvechPsi: dvec(H)/dvech(Psi)
%
% This version: August 24, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% ------------------------------------------------------

%% Compute Cholesky decomposition

H = chol(Psi, 'lower');


%% Compute derivative if desired

if nargout > 1
    
    d = size(Psi, 1); % Dimension
    
    % Index rearrangement vectors
    rearrange_vectransp = reshape(reshape(1:d^2,d,d)',d^2,1); % Vector that rearranges indices of vec(H) into indices of vec(H')
    select_vech = vec_to_vech(d); % Vector that selects linear indices of vech(H) from indices of vec(H)
    
    % Intermediate calculations, cf. Lütkepohl
    the_kron = kron(H, eye(d));
    sum_kron = the_kron + the_kron(rearrange_vectransp, :);
    sum_kron_vech = sum_kron(select_vech, :);
    
    % dvech(H)/dvech(Psi), cf. Lütkepohl
    dvechH_dvechPsi = inv(sum_kron_vech(:, select_vech));
    
    % Rearrange to obtain dvec(H)/dvech(Psi)
    dvecH_dvechPsi = zeros(d^2, size(dvechH_dvechPsi, 2));
    dvecH_dvechPsi(select_vech,:) = dvechH_dvechPsi;
    
end

end