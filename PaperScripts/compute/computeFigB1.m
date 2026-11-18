function r = computeFigB1()
%COMPUTEFIGB1 Shape of the liquid bridge during the plate pre-stretch (Appendix B).
%
%   Lubrication model of Appendix B (Spiegelberg et al. 1996): the plates
%   (radius R0 = 3 mm) separate from L0 = 1 mm until the mid-filament
%   radius equals a0 = 0.62 mm (the gap Lend, as in Fig 10). The positions
%   follow from Eqs. B6 and B9, so no simulation is needed.
%
%   Fields (lengths in m)
%     R0, L0, Lend, a0
%     zeta0      starting heights zeta0 = z/L0 of the Lagrangian points on
%                the free surface (xi0 = 1)
%     path.r, path.z   trajectories of those points, nH-by-numel(zeta0),
%                      from the cylinder (H = 1) to H = Lend/L0
%     H_profiles       plate gaps (in units of L0) at which the whole free
%                      surface is stored, from 1 to Lend/L0
%     profile.r, profile.z   free surface at those gaps, nS-by-numel(H_profiles)
%     V0, t_profiles   example plate speed (10 mm/s) and the times at which
%                      it reaches those gaps, t = (H - 1) L0/V0 [s]
%     z0curve.zeta0, z0curve.z0   z0 against zeta0 (Eq. B7)

p = defaultParameters();
geom = struct('R0', 3e-3, 'L0', 1e-3, 'zeta0', 0.5, 'xi0', 1);
geom.Lend = retractionGapForRadius(p.a0, geom);
r.R0 = geom.R0;
r.L0 = geom.L0;
r.Lend = geom.Lend;
r.a0 = p.a0;
Hend = geom.Lend/geom.L0;

% Lagrangian points on the surface, equally spaced on the initial cylinder
r.zeta0 = linspace(0, 1, 11);
H = linspace(1, Hend, 200)';
[zeta, xi] = lagrangianPosition(r.zeta0, H);
r.path.z = zeta*geom.L0;
r.path.r = xi*geom.R0;

% whole free surface at a few plate gaps
r.H_profiles = [1, 1.5, 2.5, 5, Hend];
zetaFine = linspace(0, 1, 2001);
[zeta, xi] = lagrangianPosition(zetaFine, r.H_profiles');
r.profile.z = zeta'*geom.L0;
r.profile.r = xi'*geom.R0;
r.V0 = 10e-3;   % m/s, example plate speed for the time labels
r.t_profiles = (r.H_profiles - 1)*geom.L0/r.V0;

% z0 against the starting height (Eq. B7)
r.z0curve.zeta0 = linspace(0, 1, 2001);
r.z0curve.z0 = shapeParameter(r.z0curve.zeta0);
end

function z0 = shapeParameter(zeta0)
% Eq. B7: 0 at the plates, -Inf at the midplane
z0 = (zeta0 - 1).*zeta0./(1 - 2*zeta0).^2;
end

function [zeta, xi] = lagrangianPosition(zeta0, H)
% Positions at plate gaps H (column) of the surface points that start at
% heights zeta0 (row) on the cylinder xi0 = 1: Eq. B6 for zeta, Eq. B9 for xi
z0 = shapeParameter(zeta0);
side = sign(zeta0 - 0.5);   % + above the midplane, - below
zeta = H/2.*(1 + side.*sqrt(H./(H - 4*z0)));
ratio = (1 - 4*z0./H)./(1 - 4*z0);
mid = isinf(z0);
ratio(:, mid) = repmat(1./H, 1, nnz(mid));   % limit z0 -> -Inf: xi = H^(-3/4)
xi = ratio.^(3/4);
end
