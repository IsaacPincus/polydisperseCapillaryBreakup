%CREATEFIG3 Example simulation and extraction of tau_EC (Fig 3).
%   Run this script from anywhere; results are also left in 'results'.

close all
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setupPaths();

results = computeFig3();
plotFig3(results);
