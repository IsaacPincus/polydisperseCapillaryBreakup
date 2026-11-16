%CREATEFIG9_FIG10 tau_EC vs initial radius (Fig 9) and plate retraction speed (Fig 10).
%   Run this script from anywhere; results are also left in 'results'.

close all
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setupPaths();

results = computeFig9_Fig10();
plotFig9(results);
plotFig10(results);
