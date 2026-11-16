%CREATEFIG2 Linear viscoelastic moduli of PS7 and PS16 (Fig 2).
%   Run this script from anywhere; results are also left in 'results'.

close all
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setupPaths();

results = computeFig2();
plotFig2(results);
fprintf('PS7:  Oh = %.4g, Ec = %.4g, beta = %.4g, c* = %.4g g/m^3\n', ...
    results.PS7.Oh, results.PS7.Ec, results.PS7.beta, results.PS7.cstar);
fprintf('PS16: Oh = %.4g, Ec = %.4g, beta = %.4g\n', ...
    results.PS16.Oh, results.PS16.Ec, results.PS16.beta);
