function R = bridgeRadius(H, geom)
%BRIDGERADIUS Neck radius of the liquid bridge at plate gap H (in units of L0).
%
%   R = bridgeRadius(H, geom) uses the liquid-bridge model of plateRetraction
%   with initial aspect ratio Lam0 = L0/R0 (geom fields R0, L0, zeta0, xi0).
%   For zeta0 = 0.5 this is R = R0 * H^(-3/(4 Lam0^2)).

Lam0 = geom.L0/geom.R0;
z0 = ((geom.zeta0-1)*geom.zeta0)/(1-2*geom.zeta0)^2;
R = (geom.xi0*((H/z0 - 4)./(H/z0 - H*4)).^(3/(4*Lam0^2)))*geom.R0;
end
