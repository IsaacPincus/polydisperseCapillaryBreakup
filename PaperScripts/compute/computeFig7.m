function r = computeFig7()
%COMPUTEFIG7 Stresses and Weissenberg numbers during thinning (Fig 7).
%
%   1000 wppm PS7 + 100 wppm PS16 with conformation-dependent drag.
%   Figs 7a,b: one mode per species, Ns = 30.
%   Figs 7c,d: blended distribution split into 30 modes, Ns = 1;
%   coilStretchMode is the mode whose Wi at ep_min just exceeds 0.5.

p = defaultParameters();
mw = loadMWDistributions(p);
cPS7 = 1000;
cPS16 = 100;

r.p = p;
r.twoMode = simulate(p, twoSpeciesModes(mw, cPS7, cPS16, p), 30);
r.multiMode = simulate(p, blendModes(mw, cPS7, cPS16, 30, p), 1);
% the sub-species that just touches the coil-stretch transition (Fig 7c,d):
% the smallest Wi >= 0.5 in the EC regime, at the strain-rate minimum
WiEC = r.multiMode.sol.Wi(:, r.multiMode.relax.ep_min_index);
WiEC(WiEC < 0.5) = Inf;
[~, r.multiMode.coilStretchMode] = min(WiEC);
end

function s = simulate(p, modes, Ns)
s.modes = modes;
s.pm = modeProperties(modes, p, Ns);
s.sol = solveFilamentThinning(p, s.pm, 'Drag', 'conformationDependent');
s.relax = simulatedRelaxationTime(s.sol.t, s.sol.a, s.sol.ep);
end
