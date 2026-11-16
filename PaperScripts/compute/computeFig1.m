function r = computeFig1()
%COMPUTEFIG1 MW distributions of PS7 and PS16 (Fig 1) and their 30-mode split (Fig A1).
%
%   Fig 1:  each species as a single mode, 1000 wppm PS7 + 500 wppm PS16.
%   Fig A1: the blended distribution split into 30 modes.

p = defaultParameters();
mw = loadMWDistributions(p);

r.cPS7 = 1000;
r.cPS16 = 500;
r.x = mw.x;
r.yPS7 = mw.yPS7;
r.yPS16 = mw.yPS16;
r.modes = twoSpeciesModes(mw, r.cPS7, r.cPS16, p);

r.A1 = blendModes(mw, r.cPS7, r.cPS16, 30, p);
end
