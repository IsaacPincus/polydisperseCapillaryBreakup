function r = computeFig8()
%COMPUTEFIG8 Which molecular weights set tau_EC in a polydisperse blend (Fig 8).
%
%   The blended distribution is split into 80 modes (Ns = 1,
%   conformation-dependent drag).
%   Fig 8a: 1000 wppm PS7 + 100 wppm PS16 - stress-, analytically- and
%           number-weighted MW distributions at the strain-rate minimum.
%   Fig 8b: weighted mean MW against PS16 concentration for 200, 500 and
%           1000 wppm PS7.

p = defaultParameters();
Nm = 80;
mw = loadMWDistributions(p);
r.p = p;

% --- Fig 8b: concentration sweep ------------------------------------------
r.cPS7 = [200, 500, 1000];
r.cPS16 = linspace(0, 300, 50);
nI = numel(r.cPS16);
nJ = numel(r.cPS7);
sweepFields = {'t_c', 'ep_min', 'tau_ep_min', 'tau_EC', ...
    'meanStressP', 'meanStressAna', 'meanNumber'};
for k = 1:numel(sweepFields)
    r.sweep.(sweepFields{k}) = zeros(nI, nJ);
end
for jj = 1:nJ
    for ii = 1:nI
        s = simulate(p, mw, r.cPS7(jj), r.cPS16(ii), Nm);
        m = weightedMeanMW(s);
        r.sweep.t_c(ii,jj) = s.relax.t_c;
        r.sweep.ep_min(ii,jj) = s.relax.ep_min;
        r.sweep.tau_ep_min(ii,jj) = s.relax.tau_ep_min;
        r.sweep.tau_EC(ii,jj) = s.relax.tau_EC;
        r.sweep.meanStressP(ii,jj) = m.stress;
        r.sweep.meanStressAna(ii,jj) = m.analytical;
        r.sweep.meanNumber(ii,jj) = m.number;
    end
end

% --- Fig 8a: single blend --------------------------------------------------
s = simulate(p, mw, 1000, 100, Nm);
r.single = s;
k = s.relax.ep_min_index;
M = s.modes.M;
n = s.modes.n;
stress_min_ana = 2*s.pm.g.*s.pm.lambda_m'.*s.pm.L.^2.*s.relax.ep_min;
stress_min = s.sol.stress(k,:);
xte = s.sol.t(k) - s.relax.t_c;

% MW distributions weighted by analytical stress, simulated stress, number
r.single.X1 = M*1e-6;
r.single.YMatrix1 = [(-stress_min_ana./trapz(M,stress_min_ana))';
                     -stress_min./trapz(M,stress_min);
                     -n./trapz(M,n)];
% inset: capillary, viscous + inertial, and modal polymer stresses
a = s.sol.a;
ep = s.sol.ep;
r.single.X2 = s.sol.t - s.relax.t_c;
r.single.YMatrix2 = [p.chi./a, 3*p.eta_s*ep + 1/8*p.rho*ep.^2.*a.^2, s.sol.stress];
r.single.X3 = [xte, xte];
end

function s = simulate(p, mw, cPS7, cPS16, Nm)
s.modes = blendModes(mw, cPS7, cPS16, Nm, p);
s.pm = modeProperties(s.modes, p, 1);
s.sol = solveFilamentThinning(p, s.pm, 'Drag', 'conformationDependent');
s.relax = simulatedRelaxationTime(s.sol.t, s.sol.a, s.sol.ep);
end

function m = weightedMeanMW(s)
% mean MW [MDa] weighted by the modal polymer stress at the strain-rate
% minimum, by its analytical estimate 2 g lambda L^2 ep, and by number
k = s.relax.ep_min_index;
xs = s.modes.M*1e-6;   % descending, hence the minus signs
n = s.modes.n;

stress_min = s.sol.stress(k,:);
ys = -stress_min./trapz(xs,stress_min);
m.stress = -trapz(xs,xs.*ys);

stress_ana = 2*s.pm.g.*s.pm.lambda_m'.*s.pm.L.^2.*s.sol.ep(k);
ys = -stress_ana./trapz(xs,stress_ana);
m.analytical = -trapz(xs,xs.*ys');

m.number = sum(n)/sum(n./xs);
end
