function modes = blendModes(mw, cPS7, cPS16, Nm, p)
%BLENDMODES Split the combined PS7 + PS16 distribution into Nm modes.
%
%   modes = blendModes(mw, cPS7, cPS16, Nm, p) mixes the two weight-
%   fraction distributions in proportion to their concentrations (wppm)
%   and splits the result with splitMWDist.
%
%   Fields: n, M (1xNm, descending M), molDensity, y (blended
%   distribution), buckets (see splitMWDist).

massPolymer = cPS16 + cPS7;
weightFractionPS16 = cPS16/massPolymer;
modes.y = mw.yPS7*(1-weightFractionPS16) + mw.yPS16*weightFractionPS16;
[modes.n, modes.M, modes.molDensity, modes.buckets] = ...
    splitMWDist(mw.x, modes.y, Nm, p.nu, massPolymer/p.vol_solvent);
end
