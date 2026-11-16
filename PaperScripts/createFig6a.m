%CREATEFIG6A Thinning of PS7 + PS16 blends (Fig 6a).
%   Run this script from anywhere; results are also left in 'results'.

close all
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setupPaths();

results = computeFig6a();
plotFig6a(results);
