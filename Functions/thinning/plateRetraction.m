function pull = plateRetraction(pm, V0, geom)
%PLATERETRACTION Polymer pre-stretch during CaBER plate separation.
%
%   pull = plateRetraction(pm, V0, geom) integrates the FENE-P modes of pm
%   (from modeProperties) while the plates separate at speed V0 [m/s]
%   from gap geom.L0 to geom.Lend. The filament strain rate follows the
%   liquid-bridge model with initial aspect ratio Lam0 = L0/R0 and
%   zeta0, xi0 as in the paper; plate gap H(t) = 1 + t V0/L0 (in units of L0).
%
%   geom fields: R0, L0, Lend, zeta0, xi0 (all SI). retractionGapForRadius
%   gives the Lend that ends the retraction at a chosen radius.
%
%   Output fields
%     tpull     duration of the plate motion [s]
%     Rpull     filament radius when the plates stop [m]
%     epend     strain rate when the plates stop [1/s]
%     conformation  state vector at the end (see packConformation), the
%                   'InitialConformation' for solveFilamentThinning
%     stress    polymer stress difference of each species at the end
%               [Pa], 1-by-Nm (g f sum_s (Azz - Arr))

[Nm, Ns] = size(pm.lambda_ms);
Lam0 = geom.L0/geom.R0;
Hend = geom.Lend/geom.L0;
Hfunc = @(t) 1 + t;
Hderivfunc = @(t) 1;
z0 = ((geom.zeta0-1)*geom.zeta0)/(1-2*geom.zeta0)^2;

pull.tpull = fzero(@(t) Hfunc(t) - Hend, 1)*geom.L0/V0;

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
y0 = ones(Nm*Ns*2, 1);
sol = ode15s(@(t,y) pullRHS(t, y, pm, Hfunc, Hderivfunc, z0, Lam0, V0, geom.L0), ...
    [0, pull.tpull], y0, options);
points = 1000;
y = deval(sol, linspace(0, pull.tpull, points));
pull.conformation = y(:,end);

H = Hfunc(pull.tpull*V0/geom.L0);
dH = Hderivfunc(pull.tpull*V0/geom.L0);
pull.Rpull = bridgeRadius(H, geom);
pull.epend = (-3/(2*Lam0^2)*1./(H/(4*z0) - 1).*dH./H)*V0/geom.L0;

[Az, Ar] = unpackConformation(pull.conformation, Nm, Ns);
L2 = (pm.L).^2;
f = L2./(L2+3*Ns-sum(Az+2*Ar,2));   % Eq. 22
pull.stress = (pm.g.*f.*sum(Az-Ar,2))';
end

function dydt = pullRHS(t, y, pm, Hfunc, Hderivfunc, z0, Lam0, V0, L0)
H = Hfunc(t*V0/L0);
dH = Hderivfunc(t*V0/L0);
ep = (-3/(2*Lam0^2)*1./(H/(4*z0) - 1).*dH./H)*V0/L0;

[Nm, Ns] = size(pm.lambda_ms);
[Az, Ar] = unpackConformation(y, Nm, Ns);
L2 = (pm.L).^2;
f = L2./(L2+3*Ns-sum(Az+2*Ar,2));   % Eq. 22
Azdot = 2*ep*Az - f./pm.lambda_ms.*(Az-1);
Ardot = -ep*Ar - f./pm.lambda_ms.*(Ar-1);
dydt = packConformation(Azdot, Ardot);
end
