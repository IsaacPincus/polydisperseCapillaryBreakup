function Lend = retractionGapForRadius(a, geom)
%RETRACTIONGAPFORRADIUS Final plate gap at which the filament radius is a.
%
%   Lend = retractionGapForRadius(a, geom) inverts bridgeRadius: the plates
%   must stop at gap Lend [m] for the filament to start thinning from
%   radius a [m]. The result does not depend on the retraction speed.

Lend = fzero(@(L) bridgeRadius(L/geom.L0, geom) - a, geom.L0*[1, 10]);
end
