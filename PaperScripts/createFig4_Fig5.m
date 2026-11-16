%CREATEFIG4_FIG5 Single-species thinning vs experiment (Fig 4) and tau_EC/tau_Z vs c/c* (Fig 5).
%   Run this script from anywhere; results are also left in 'results'.

close all
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setupPaths();

results = computeFig4_Fig5();
plotFig4(results);
plotFig5(results);
