%CREATEFIG6B tau_EC of PS7 + PS16 blends vs PS16 concentration (Fig 6b).
%   Run this script from anywhere; results are also left in 'results'.

close all
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setupPaths();

results = computeFig6b();
plotFig6b(results);
