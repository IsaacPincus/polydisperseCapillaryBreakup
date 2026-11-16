function out = simulateConcentrationSweep(p, mw, cPS7, cPS16, Ns, drag)
%SIMULATECONCENTRATIONSWEEP Thinning simulations for a list of PS7 + PS16 blends.
%
%   out = simulateConcentrationSweep(p, mw, cPS7, cPS16, Ns, drag) runs
%   one simulation per entry of cPS7/cPS16 (wppm), each species being a
%   single MW mode with Ns Zimm sub-modes, and extracts tau_EC with
%   simulatedRelaxationTime (unshifted fit).
%
%   This is the 'get_outputs' function of the original Fig 4-6 scripts.
%   drag is passed to solveFilamentThinning ('constant' or
%   'conformationDependent').
%
%   Output fields (1-by-nRuns or 1-by-nRuns cells): t_c, tau_ep_min,
%   tau_EC, aModel, epModel, tModel, stressModel, f_min (nRuns-by-2: the
%   FENE-PM factor f of [PS16, PS7] at the strain-rate minimum), and
%   lambda_m (1x2: [PS16, PS7]).

nRuns = numel(cPS16);
out.t_c = zeros(1, nRuns);
out.tau_ep_min = zeros(1, nRuns);
out.tau_EC = zeros(1, nRuns);
out.aModel = cell(1, nRuns);
out.epModel = cell(1, nRuns);
out.tModel = cell(1, nRuns);
out.stressModel = cell(1, nRuns);
out.f_min = zeros(nRuns, 2);

for ii = 1:nRuns
    modes = twoSpeciesModes(mw, cPS7(ii), cPS16(ii), p);
    pm = modeProperties(modes, p, Ns);
    sol = solveFilamentThinning(p, pm, 'Drag', drag);
    r = simulatedRelaxationTime(sol.t, sol.a, sol.ep);

    out.t_c(ii) = r.t_c;
    out.tau_ep_min(ii) = r.tau_ep_min;
    out.tau_EC(ii) = r.tau_EC;
    out.aModel{ii} = sol.a;
    out.epModel{ii} = sol.ep;
    out.tModel{ii} = sol.t;
    out.stressModel{ii} = sol.stress;
    out.f_min(ii,:) = sol.f(r.ep_min_index,:);
end
out.lambda_m = pm.lambda_m;
end
