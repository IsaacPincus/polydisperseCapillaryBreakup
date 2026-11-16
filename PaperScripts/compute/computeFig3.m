function r = computeFig3()
%COMPUTEFIG3 Example simulation and how tau_EC is extracted (Fig 3).
%
%   1000 wppm PS7 + 500 wppm PS16, one mode per species with 30 Zimm
%   sub-modes, constant drag. tau_EC from an exponential fit of a(t - t_c)
%   where ep is within 10% of ep_min.

p = defaultParameters();
Ns = 30;
mw = loadMWDistributions(p);

r.modes = twoSpeciesModes(mw, 1000, 500, p);
r.pm = modeProperties(r.modes, p, Ns);
r.sol = solveFilamentThinning(p, r.pm);
r.relax = simulatedRelaxationTime(r.sol.t, r.sol.a, r.sol.ep);
end
