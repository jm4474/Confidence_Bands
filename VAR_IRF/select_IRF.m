function sel = select_IRF(d, nshock, hori, vari, shock)
% ------------------------------------------------
% Selection vector "sel" for single IRF
% such that IRF of interest is given by Theta(sel)
%
% Inputs:
% - d:      number of variables
% - nshock: number of identified shocks
% - hori:   number of horizons in Theta
% - vari:   index of response variable to select
% - shock:  index of shock to select
% Outputs:
% - sel:    selection vector
%
% This version: August 24, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% ------------------------------------------------


%% Selection vector

irf_zeros = zeros(d, nshock*(hori+1));
irf_zeros(vari, nshock*(0:hori)+shock) = 1;
sel = (irf_zeros(:)==1);


end