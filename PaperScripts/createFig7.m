%CREATEFIG7 Stresses and effective Weissenberg numbers during thinning (Fig 7).
%   Run this script from anywhere; results are also left in 'results'.

close all
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setupPaths();

results = computeFig7();
plotFig7(results);
