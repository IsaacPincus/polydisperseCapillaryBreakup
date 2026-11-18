function rhoEff = inertialDensity(rho, alpha_C1)
%INERTIALDENSITY Fluid density scaled for the 0D model's inertial term.
%
%   rhoEff = inertialDensity(rho, alpha_C1) returns rho*C1*2 with
%   C1 = 9/(4*alpha_C1^3), so that the solver's inertial term
%   rhoEff a^2 ep^2/8 is (C1/4) rho a^2 ep^2 as in Eq. 8 (and Zinelis et
%   al. 2024). The original scripts used rho*C1*8, i.e. 4 times more.

C1 = 9/(4*alpha_C1^3);
rhoEff = rho*C1*2;
end
