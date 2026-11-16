function M = readDataFile(fileName, sheet, columns)
%READDATAFILE Read a numeric matrix from a file in the Data folder.
%
%   M = readDataFile(fileName)                 whole file
%   M = readDataFile(fileName, sheet)          one spreadsheet sheet
%   M = readDataFile(fileName, sheet, columns) only the listed columns
%
%   Uses detectImportOptions/readmatrix exactly as the original scripts
%   did, so header handling is unchanged.

file = fullfile(dataDirectory(), fileName);
if nargin > 1 && ~isempty(sheet)
    opts = detectImportOptions(file, 'Sheet', sheet);
else
    opts = detectImportOptions(file);
end
if nargin > 2 && ~isempty(columns)
    opts.SelectedVariableNames = columns;
end
M = readmatrix(file, opts);
end
