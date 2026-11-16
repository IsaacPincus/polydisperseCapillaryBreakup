function [t, D] = readThinningCurves(sheet, nCurves, skipFirstRow)
%READTHINNINGCURVES Experimental CaBER curves from the thinning spreadsheet.
%
%   [t, D] = readThinningCurves(sheet, nCurves) reads the first nCurves
%   (t - t_c [s], D [mm]) column pairs of the given sheet of
%   ThinningCurves_PS7_PS16_andMixtures.xlsx. Pairs sit every 4 columns.
%   Each column of t and D is one curve (NaN-padded).
%
%   skipFirstRow (default true) drops the first data row, as most of the
%   original scripts did. Fig 4b read the PS7-only curves without
%   dropping it, so it passes false.

if nargin < 3
    skipFirstRow = true;
end

columns = reshape([1:4:4*nCurves; 2:4:4*nCurves], 1, []);
M = readDataFile('ThinningCurves_PS7_PS16_andMixtures.xlsx', sheet, columns);
if skipFirstRow
    M = M(2:end, :);
end
t = M(:, 1:2:end);
D = M(:, 2:2:end);
end
