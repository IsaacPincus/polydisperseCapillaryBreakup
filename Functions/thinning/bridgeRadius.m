function R = bridgeRadius(H, geom)
%BRIDGERADIUS Neck radius of the liquid bridge at plate gap H (in units of L0).
%
%   R = bridgeRadius(H, geom) is Eq. B9 of the lubrication model used by
%   plateRetraction (geom fields R0, zeta0, xi0). It does not depend on the
%   aspect ratio L0/R0. For zeta0 = 0.5 this is R = R0 * H^(-3/4).

z0 = ((geom.zeta0-1)*geom.zeta0)/(1-2*geom.zeta0)^2;
R = (geom.xi0*((H/z0 - 4)./(H/z0 - H*4)).^(3/4))*geom.R0;
end
