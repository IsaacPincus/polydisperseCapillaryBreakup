function Ub = ubcal(hstar, nkperb)
%UBCAL Drag of an isolated Zimm blob relative to a Stokes sphere.
%
%   Written by R. Prabhakar (Prabhakar Ranganathan), who provided the
%   conformation-dependent drag implementation used in the paper. See
%   R. Prabhakar et al., J. Rheol. 60, 345-366 (2016), and
%   R. Prabhakar et al., Phys. Rev. Fluids 2, 011301 (2017).
%
%   Ub = ubcal(hstar, nkperb) returns U_b = zeta_b/(3 pi eta_s xi_b) for a
%   blob of nkperb Kuhn segments with hydrodynamic interaction parameter
%   hstar (see the discussion on p. 350 of the J. Rheol. paper, left column).
%   Works elementwise on nkperb.

aK = pi^0.5*hstar;                           % hydrodynamic radius of a Kuhn segment
X = nkperb.^0.5;
Pr = 0.9119 - 0.0191./X + 0.1009./X.^2;      % fit of R/N^0.5 vs 1/N^0.5
Psr = 1.6366 - 1.1434./X - 0.4521./X.^2;     % fit of S R/N^2 vs 1/N^0.5 (KR double sum of 1/r_ij)
Ub = 2*aK./(Pr./X + 2*aK*Psr);
end
