%CREATEFIG8 MW weighting of tau_EC in polydisperse blends (Fig 8). Takes the longest (~20 min).
%   Run this script from anywhere; results are also left in 'results'.

close all
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setupPaths();

results = computeFig8();
plotFig8(results);
