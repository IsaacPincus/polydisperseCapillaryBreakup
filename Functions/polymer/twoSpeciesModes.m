function modes = twoSpeciesModes(mw, cPS7, cPS16, p)
%TWOSPECIESMODES One mode per polymer species for a PS7 + PS16 blend.
%
%   modes = twoSpeciesModes(mw, cPS7, cPS16, p) represents each species by
%   a single mode (splitMWDist with N = 1) and combines them by number.
%   Mode 1 is PS16, mode 2 is PS7. Concentrations are in wppm.
%
%   Fields: n (1x2 number fractions), M (1x2 MW), molDensity (mol/m^3),
%   and the per-species values nPS7, MPS7, molDensityPS7 (same for PS16).

[modes.nPS7, modes.MPS7, modes.molDensityPS7] = ...
    splitMWDist(mw.x, mw.yPS7, 1, p.nu, cPS7/p.vol_solvent);
[modes.nPS16, modes.MPS16, modes.molDensityPS16] = ...
    splitMWDist(mw.x, mw.yPS16, 1, p.nu, cPS16/p.vol_solvent);

modes.molDensity = modes.molDensityPS7 + modes.molDensityPS16;
molFracPS7 = modes.molDensityPS7/modes.molDensity;
modes.n = [modes.nPS16*(1-molFracPS7), modes.nPS7*molFracPS7];
modes.M = [modes.MPS16, modes.MPS7];
end
