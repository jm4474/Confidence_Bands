function Phi = MARep(AL,hori)
% -------------------------------------------------------------------------
% Transforms the A(L) parameters of a reduced-form VAR
% into the coefficients Phi of the MA representation.
% 
% Inputs:
% - AL:     VAR model coefficients
% - hori:   largest response horizon
% Outputs:
% - Phi:    MA representation coefficients
% 
% This version: August 24, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% -------------------------------------------------------------------------


%% Reshape AL into a 3-D array

d = size(AL,1);
tau = size(AL,2)/d;
vecAL = reshape(AL,[d,d,tau]); 


%% Initialize the value of the auxiliary array vecALrevT

vecALrevT = zeros(d,d,hori);
for i=1:hori
    if i<(hori-tau)+1
        vecALrevT(:,:,i) = zeros(d,d);
    else
        vecALrevT(:,:,i) = vecAL(:,:,(hori-i)+1)';
    end
end
vecALrevT = reshape(vecALrevT,[d,d*hori]);


%% MA coefficients

Phi = repmat(vecAL(:,:,1),[1,hori]);
for i=1:hori-1
    Phi(:,(d*i)+1:(d*(i+1))) = [eye(d),Phi(:,1:d*i)] * vecALrevT(:,(hori*d-(d*(i+1)))+1:end)';  
end

end

