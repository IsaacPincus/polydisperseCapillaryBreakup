function ax = tauEcFigure(r, x, tau_EC, xLabel)
%TAUECFIGURE tau_EC against x, one line with markers per PS16 concentration,
%   with dotted lines at the Zimm times of PS16 and PS7 (Figs 9 and 10a).
%   The y-limits are the same for both figures. No legend.

paperFigure();
ax = gca;
colours = winter(numel(r.cPS16));
for k = 1:numel(r.cPS16)
    plot(x, tau_EC(:,k), '-o', 'Color', colours(k,:), 'MarkerFaceColor', colours(k,:), ...
        'DisplayName', sprintf('$c_\\mathrm{PS16} = %d$ [wppm]', r.cPS16(k)));
end
yline(r.tauZ(1), 'k:', '$\tau_{Z,\mathrm{PS16}}$', 'Interpreter', 'latex', ...
    'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', ...
    'FontSize',16, 'LineWidth', 2, 'Color',[11, 140, 14]/255);
yline(r.tauZ(2), 'k:', '$\tau_{Z,\mathrm{PS7}}$', 'Interpreter', 'latex', ...
    'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', ...
    'FontSize',16, 'LineWidth', 2, 'Color','b');
set(ax, 'XScale', 'log', 'XMinorTick', 'on')
ylim(ax, [0, 1.1*max([r.radius.tau_EC(:); r.retraction.tau_EC(:); r.tauZ(:)])])
xlabel(xLabel)
ylabel('$\tau_\mathrm{EC}$ [s]')
end
