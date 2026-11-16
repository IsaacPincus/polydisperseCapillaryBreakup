function rhoEff = inertialDensity(rho, alpha_C1)
%INERTIALDENSITY Fluid density scaled for the 0D model's inertial term.
%
%   rhoEff = inertialDensity(rho, alpha_C1) returns rho*C1*8 with
%   C1 = 9/(4*alpha_C1^3), the form used throughout the paper scripts.

C1 = 9/(4*alpha_C1^3);
rhoEff = rho*C1*8;
end
