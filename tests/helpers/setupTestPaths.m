function setupTestPaths()
%SETUPTESTPATHS Put the refactored project and the test helpers on the path.
testsDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fileparts(testsDir));
setupPaths();
addpath(fullfile(testsDir, 'helpers'));
end
