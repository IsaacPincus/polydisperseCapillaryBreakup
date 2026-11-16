function mw = loadMWDistributions(p)
%LOADMWDISTRIBUTIONS PS7 and PS16 weight-fraction distributions on one grid.
%
%   mw = loadMWDistributions(p) reads PS7MWDist.csv and PS16MWDist.csv and
%   interpolates both onto a common, uniformly spaced MW grid with
%   p.noPoints points (method p.interpMethod, zero outside the data).
%
%   Fields: x (MW grid, g/mol), yPS7, yPS16 (weight fraction per g/mol),
%   and the raw data PS7MW, PS7NormMass, PS16MW, PS16NormMass.

dataPS7 = readDataFile('PS7MWDist.csv');
dataPS16 = readDataFile('PS16MWDist.csv');

mw.PS7MW = dataPS7(:,1);
mw.PS7NormMass = dataPS7(:,2);
mw.PS16MW = dataPS16(:,1);
mw.PS16NormMass = dataPS16(:,2);

allPS_data = [mw.PS7MW; mw.PS16MW];
mw.x = linspace(min(allPS_data), max(allPS_data), p.noPoints);
mw.yPS7 = interp1(mw.PS7MW, mw.PS7NormMass, mw.x, p.interpMethod, 0);
mw.yPS16 = interp1(mw.PS16MW, mw.PS16NormMass, mw.x, p.interpMethod, 0);
end
