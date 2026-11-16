%CREATEFIG1_FIGA1 MW distributions: Fig 1 and appendix Fig A1.
%   Run this script from anywhere; results are also left in 'results'.

close all
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setupPaths();

results = computeFig1();
plotFig1(results);
