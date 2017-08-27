function PosteriorDraws = UhligPosterior(vecnuAL,vecPsi,X,numdraws)
% ------------------------------------------------------
% Gnerates draws from the Normal-Wishart posterior
% specified in Uhlig (2005)
%
% Inputs:
% - vecnuAL:        vectorized [nuhat,ALhat];
% - vecPsi:         vectorized Psihat
% - X:              matrix of VAR covariates (T x (1 + d*tau)) (includes the constant term)
% Outputs:
% - PosteriorDraws: struct with following fields
% PosteriorDraws fields:
% - A:              draws of coefficient matrices
% - Psi:            draws of innovation variance matrices
%
% This version: August 24, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% ------------------------------------------------------


%% 0) Define values for the simulation

    d = sqrt(length(vecPsi));
    tau = (size(X,2)-1)/d;
    T = size(X,1);
    S_X = (X'*X)/T;
    

%% 1) Take I draws from the random variables Z and W 

    Z=randn(d,1,T,numdraws);
    W=randn(d*(size(X,2)),1,numdraws);

    
%% 2) Generate draws of Psi conditional on the data    

    Aux1=mean(bsxfun(@times,Z,permute(Z,[2,1,3,4])),3); %This is [n,n,1,Draws];
    Psi=reshape(vecPsi,[d,d]);
    PsiP=zeros(d,d,numdraws);
    for ix=1:numdraws
        PsiP(:,:,ix)=(Psi^(1/2))*(Aux1(:,:,1,ix)\(Psi^(1/2)));
    end
    
    
%% 3) Generate draws from vecnuA conditional on Psi and the data
   
   vecnuApaux=zeros(d*(size(X,2)),numdraws);
   for ix=1:numdraws
       vecnuApaux(:,ix)=vecnuAL + chol((kron(S_X\eye(size(X,2)),PsiP(:,:,ix)/T)))'*W(:,ix);
   end
   
   vecAp=reshape(vecnuApaux(d+1:end,:),[d,d*tau,numdraws]);
   
   PosteriorDraws.A=vecAp;
   PosteriorDraws.Psi=PsiP;

end

