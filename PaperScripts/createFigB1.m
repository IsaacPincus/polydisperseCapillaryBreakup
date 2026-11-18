%CREATEFIGB1 Liquid-bridge shape during the plate pre-stretch (Appendix B).
%   Run this script from anywhere; results are also left in 'results'.

close all
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setupPaths();

results = computeFigB1();
plotFigB1(results);
