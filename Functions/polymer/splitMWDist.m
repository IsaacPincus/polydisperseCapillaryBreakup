function [n, M, molDensity, buckets] = splitMWDist(x, y, N, nu, massDensity)
%SPLITMWDIST Discretise a MW distribution into N modes.
%
%   [n, M, molDensity, buckets] = splitMWDist(x, y, N, nu, massDensity)
%
%   Splits the weight-fraction distribution y(x) into N buckets of equal
%   width on the grid x. For each bucket it returns the extensibility-
%   averaged MW, M = (int y x^(1+nu) / int y)^(1/(1+nu)), and its number
%   fraction n (normalised to sum to 1). Modes are sorted by descending M.
%
%   Inputs
%     x, y         MW grid [g/mol] and weight fraction per g/mol
%     N            number of modes
%     nu           exponent in the extensibility average
%     massDensity  polymer mass per volume of solvent [g/m^3]
%
%   Outputs
%     n, M         1-by-N number fractions and molecular weights
%     molDensity   total polymer number density [mol/m^3]
%     buckets      unsorted per-bucket data for plotMWDistSplit:
%                  edges (grid indices), M, n

noPoints = length(x);
molDensity = trapz(x, y*massDensity./x);

edges = round(linspace(1, noPoints, N+1));
ML = zeros(1, N);
ni = zeros(1, N);
for ii = 1:N
    range = edges(ii):edges(ii+1);
    yii = zeros(size(y));
    yii(range) = y(range);
    ML(ii) = (trapz(x, yii.*x.^(1+nu))/trapz(x, yii))^(1/(1+nu));
    if isnan(ML(ii))
        ML(ii) = mean(x(range));
    end
    ni(ii) = trapz(x, yii./x);
end
total_n = trapz(x, y./x);
n = ni/total_n;

buckets.edges = edges;
buckets.M = ML;
buckets.n = n;

[M, order] = sort(ML, 'descend');
n = n(order);
end
