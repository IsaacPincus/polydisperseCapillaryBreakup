function root = setupPaths()
%SETUPPATHS Add the project's code folders to the MATLAB path.
%
%   root = setupPaths() adds Functions/ and PaperScripts/ (with
%   subfolders) and returns the project root folder. Data files are read
%   by absolute path, so the working folder does not matter.

root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, 'Functions')));
addpath(genpath(fullfile(root, 'PaperScripts')));
% ode15s warns about near-singular Jacobians during the final collapse
warning('off', 'MATLAB:nearlySingularMatrix');
end
