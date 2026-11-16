function r = computeFig9_Fig10()
%COMPUTEFIG9_FIG10 tau_EC vs initial radius (Fig 9) and vs plate retraction
%   speed, with the polymer stress at the end of the retraction (Fig 10).
%
%   1000 wppm PS7 with 1-1000 wppm PS16, one mode per species, Ns = 10,
%   constant drag.
%   radius:     a0 from 0.1 to 10 mm, filament starts at equilibrium.
%   retraction: plates pulled apart at V0 = 0.1-2 mm/s from L0 = 1 mm to
%               the gap Lend at which the filament radius equals p.a0
%               (0.62 mm, as in the other figures); thinning starts from
%               the pre-stretched state. stress16, stress7 are the polymer
%               stresses of PS16 and PS7 when the plates stop [Pa].
%   tauZ holds the Zimm times [PS16, PS7] marked in the figures; a0 is the
%   radius at the end of the retraction.

p = defaultParameters();
a0 = p.a0;   % the sweeps below overwrite p.a0
r.a0 = a0;
Ns = 10;
mw = loadMWDistributions(p);
r.cPS7 = 1000;
r.cPS16 = [1, 50, 100, 200, 500, 1000];
cPS7 = r.cPS7;

% --- changing the initial radius (Fig 9) ----------------------------------
r.radius.a0_vals = logspace(-1, 1, 10)*1e-3;
for jj = 1:numel(r.cPS16)
    pm = modeProperties(twoSpeciesModes(mw, cPS7, r.cPS16(jj), p), p, Ns);
    for ii = 1:numel(r.radius.a0_vals)
        p.a0 = r.radius.a0_vals(ii);
        sol = solveFilamentThinning(p, pm);
        relax = simulatedRelaxationTime(sol.t, sol.a, sol.ep);
        r.radius.tau_ep_min(ii,jj) = relax.tau_ep_min;
        r.radius.tau_EC(ii,jj) = relax.tau_EC;
    end
end

% --- changing the plate retraction speed (Fig 10) -------------------------
geom.R0 = 3e-3;     % plate radius
geom.L0 = 1e-3;     % initial gap
geom.zeta0 = 0.5;
geom.xi0 = 1;
geom.Lend = retractionGapForRadius(a0, geom);   % final gap: radius a0
r.retraction.Lend = geom.Lend;
r.retraction.V0_vals = linspace(0.1, 2, 10)*1e-3;
for jj = 1:numel(r.cPS16)
    pm = modeProperties(twoSpeciesModes(mw, cPS7, r.cPS16(jj), p), p, Ns);
    for ii = 1:numel(r.retraction.V0_vals)
        pull = plateRetraction(pm, r.retraction.V0_vals(ii), geom);
        p.a0 = pull.Rpull;
        sol = solveFilamentThinning(p, pm, ...
            'InitialConformation', pull.conformation);
        relax = simulatedRelaxationTime(sol.t, sol.a, sol.ep);
        r.retraction.Rpull(ii,jj) = pull.Rpull;
        r.retraction.epend(ii,jj) = pull.epend;
        r.retraction.stress16(ii,jj) = pull.stress(1);
        r.retraction.stress7(ii,jj) = pull.stress(2);
        r.retraction.tau_ep_min(ii,jj) = relax.tau_ep_min;
        r.retraction.tau_EC(ii,jj) = relax.tau_EC;
    end
end
r.retraction.lastRun = sol;
r.tauZ = pm.lambda_m;   % same for every blend
end
