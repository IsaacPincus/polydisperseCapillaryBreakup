function plotFig3(r)
%PLOTFIG3 Exponential fit used for tau_EC (Fig 3).
%   Draws the radius with the fit and the strain rate as separate figures,
%   then Fig 3 as in the paper: the radius plot with the strain rate as an
%   inset (top right), with arrows to t_c at max(ep) and to ep_min.
%   Axis limits follow the simulation: time runs from just before t_c to
%   breakup; dashed blue ticks mark the ends of the fitted EC window. Short
%   horizontal lines at the left of the radius plot mark the initial radius
%   a0 and the fit prefactor a1 (a = a1 exp(-(t - t_c)/(3 tau_EC))).

t = r.sol.t;
a = r.sol.a;
ep = r.sol.ep;
rx = r.relax;
tm = t - rx.t_c;
tLimits = [-0.05, tm(end)];
fitEnds = rx.indices_EC([1, end]);

radiusLegend = {'Exponential Fit', 'Simulation', 'Points used for fit'};

% --- radius and fit -------------------------------------------------------
figure();
axRadius = gca;
hold on
set(gca, 'YScale', 'log')
fitPlot = rx.fitobject;
fitPlot.a = fitPlot.a*1e3;      % m -> mm
plot(fitPlot)
plot(tm, a*1e3, 'k');
plot(tm(rx.indices_EC), a(rx.indices_EC)*1e3, 'or');
legend(radiusLegend)
xlabel('$t-t_c$ [s]');
ylabel('$a$ [mm]');
xlim(tLimits);
ylim([1e-3, 1]);
fitWindowTicks(tm(fitEnds), a(fitEnds)*1e3, diff(tLimits));
radiusMark(a(1)*1e3, '$a_0$', tLimits);
radiusMark(rx.fitobject.a*1e3, '$a_1$', tLimits);

% --- strain rate ----------------------------------------------------------
figure();
axRate = gca;
hold on
set(gca, 'YScale', 'log')
plot(tm, ep, 'k');
plot(rx.t_ep_min - rx.t_c, rx.ep_min, 'bd');
plot(tm, ones(size(ep))*rx.ep_min, 'b--', 'linewidth', 0.5);
plot(tm, ones(size(ep))*rx.ep_min*rx.threshold, 'b--', 'linewidth', 0.5);
legend({'Simulation', 'Min $\dot{\epsilon}$'})
xlabel('$t-t_c$ [s]');
ylabel('$\dot{\epsilon}$ [s$^{-1}$]');
xlim(tLimits);
ylim([rx.ep_min/2, 100*rx.ep_min]);

% --- Fig 3: radius with the strain rate as an inset ----------------------
fig = figure();
ax = copyobj(axRadius, fig);   % legends are not copied with the axes
legend(ax, radiusLegend, 'Location', 'southwest')
inset = copyobj(axRate, fig);
inset.Position = [0.55 0.55 0.33 0.32];
% decade y-limits so the inset has labelled ticks at each power of 10
inset.YLim = 10.^[floor(log10(rx.ep_min/2)), ceil(log10(max(ep)))];
inset.YTick = 10.^(log10(inset.YLim(1)):log10(inset.YLim(2)));
inset.YMinorTick = 'on';
drawnow   % final inset position and limits, for the arrows
ep_c = ep(find(tm >= 0, 1));   % strain rate at t_c (its first maximum)
span = diff(tLimits);
insetArrow(inset, [0.25*span, 0], [ep_c, ep_c], ...
    '$t_c$ at $\dot{\epsilon}_\mathrm{max}$', 'left', 'middle');
x_min = rx.t_ep_min - rx.t_c;
insetArrow(inset, [x_min, x_min], [8*rx.ep_min, 1.15*rx.ep_min], ...
    '$\dot{\epsilon}_\mathrm{min}$', 'center', 'bottom');
end

function insetArrow(ax, x, y, label, hAlign, vAlign)
% text arrow from (x(1), y(1)) to (x(2), y(2)) in the data units of axes ax
% (linear x, log y); the label sits at the tail
pos = ax.Position;
xl = ax.XLim;
yl = log10(ax.YLim);
xFig = pos(1) + pos(3)*(x - xl(1))/diff(xl);
yFig = pos(2) + pos(4)*(log10(y) - yl(1))/diff(yl);
annotation(ax.Parent, 'textarrow', xFig, yFig, 'String', label, ...
    'HorizontalAlignment', hAlign, 'VerticalAlignment', vAlign, 'FontSize', 16);
end

function radiusMark(aLevel, label, tLimits)
% short horizontal line at a = aLevel from the left edge of the axes,
% labelled just after its end
x = tLimits(1) + [0, 0.3]*diff(tLimits);
plot(x, aLevel*[1, 1], 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
text(x(2) + 0.01*diff(tLimits), aLevel, label, 'VerticalAlignment', 'middle',...
    'FontSize', 14);
end

function fitWindowTicks(t, a, tSpan)
% short dashed diagonal ticks through the first and last fitted points
dt = 0.02*tSpan;
for k = 1:numel(t)
    plot(t(k) + [-dt, dt], a(k)*[10^-0.25, 10^0.25], 'b--', 'LineWidth', 0.5, ...
        'HandleVisibility', 'off');
end
end
