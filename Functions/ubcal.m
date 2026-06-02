function [Ub] = ubcal(hstar,nkperb)
%% 
% Calculates Ub, the ratio of an isotropic Zimm blob's friction coefficient 
% to the Stokes drag of a sphere of the same size.
% See discussion on page 350 of JoR paper, left column. 
% Ub = \alpha/ (3 \pi), where \alpha = \zeta/ (\eta_s \xi); 

%%
aK = pi^0.5*hstar; % hydrodynamic radius of Kuhn segment
X = nkperb.^0.5;
Pr = 0.9119 - 0.0191./X + 0.1009./X.^2; % Polynomial fit of R/N^0.5 vs 1/N^0.5
Psr = 1.6366 - 1.1434./X -0.4521./X.^2; % Polynomial fit of SR/N^2 vs 1/N^0.5 where S is the KR double sum of 1/r_ij
Ub = 2*aK./(Pr./X + 2*aK*Psr); % U_b = zeta_b/(3 pi etas xi_b) 
%Ub = 2*0.3055;

end