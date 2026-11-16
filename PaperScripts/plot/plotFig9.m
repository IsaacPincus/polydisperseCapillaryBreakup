function plotFig9(r)
%PLOTFIG9 tau_EC vs initial radius a0 (Fig 9), on the same y-axis limits
%   as Fig 10a. Dotted lines mark the Zimm times of PS16 and PS7. The
%   region Bo > 1 is shaded light grey, with Bo = rho g a0 L/chi and L the
%   gap at the end of the plate retraction. The a0 for Bo = 1 is printed,
%   together with the one for Bo = rho g a0^2/chi. The legend is in Fig 10a.

tauEcFigure(r, r.radius.a0_vals*1000, r.radius.tau_EC, '$a_0$ [mm]');
p = defaultParameters();
gravity = 9.81;   % m/s^2
a0_BoL = p.chi/(p.rhoFluid*gravity*r.retraction.Lend);   % rho g a0 L/chi = 1
a0_Bo2 = sqrt(p.chi/(p.rhoFluid*gravity));               % rho g a0^2/chi = 1
fprintf('Bo = 1 at a0 = %.3g mm (rho g a0 L/chi, L = %.3g mm) and %.3g mm (rho g a0^2/chi)\n', ...
    a0_BoL*1e3, r.retraction.Lend*1e3, a0_Bo2*1e3);
bondRegion(gca, a0_BoL*1000);
lg = legend('Location', 'southeast', 'NumColumns', 2);
w = 0.4;
annotation('textbox', [sum(lg.Position([1 3]))-w, sum(lg.Position([2 4]))+0.01, w, 0.07], ...
    'String', sprintf('$c_\\mathrm{PS7} = %g$ [wppm]', r.cPS7), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
end

function bondRegion(ax, xBo)
% light grey patch over x > xBo (Bo > 1), behind the curves, with 'Bo > 1'
% written vertically along its left edge
xline(xBo, 'k:', ' ', 'Interpreter', 'latex', ...
    'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', ...
    'FontSize',16, 'LineWidth', 0.5, 'Color',[0.1,0.1,0.1]);
xl = ax.XLim;
yl = ax.YLim;
hp = patch(ax, [xBo, xl(2), xl(2), xBo], yl([1 1 2 2]), 0.93*[1 1 1], ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
uistack(hp, 'bottom');
set(ax, 'Layer', 'top', 'XLim', xl, 'YLim', yl)   % keep ticks visible and limits fixed
text(ax, 1.03*xBo, mean(yl)*0.9, '$\mathrm{Bo} > 1$', 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top',...
    'FontSize', 14);
end
