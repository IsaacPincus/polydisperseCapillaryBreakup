function plotThinningCurves(t_exp, y_exp, sim, labels, colours, yScale)
%PLOTTHINNINGCURVES Measured (symbols) and simulated (lines) thinning curves.
%
%   t_exp, y_exp  cell arrays of experimental times [s] and radius or
%                 diameter [mm], one per curve
%   sim           output of simulateConcentrationSweep
%   labels        legend entry per curve
%   colours       one row per curve
%   yScale        multiplies the simulated radius [m] to match y_exp
%                 (1e3 for radius in mm, 2e3 for diameter in mm)
paperFigure();
for ii = 1:numel(t_exp)
    plot(t_exp{ii}, y_exp{ii}, 's:', 'Color', colours(ii,:), 'DisplayName', labels{ii})
    plot(sim.tModel{ii} - sim.t_c(ii), sim.aModel{ii}*yScale, '-', 'linewidth', 2, ...
        'Color', colours(ii,:), 'HandleVisibility', 'off');
end
set(gca, 'YScale', 'log')
xlabel('$t-t_c$ [s]');
end
