function plotFig3(r)
%PLOTFIG3 Exponential fit used for tau_EC (Fig 3a) and strain rate (Fig 3b).
%   Axis limits follow the simulation: time runs from just before t_c to
%   breakup; dashed blue ticks mark the ends of the fitted EC window.

t = r.sol.t;
a = r.sol.a;
ep = r.sol.ep;
rx = r.relax;
tm = t - rx.t_c;
tLimits = [-0.015, tm(end)];
fitEnds = rx.indices_EC([1, end]);

% --- Fig 3a: radius and fit -----------------------------------------------
figure();
hold on
set(gca, 'YScale', 'log')
fitPlot = rx.fitobject;
fitPlot.a = fitPlot.a*1e3;      % m -> mm
plot(fitPlot)
plot(tm, a*1e3, 'k');
plot(tm(rx.indices_EC), a(rx.indices_EC)*1e3, 'or');
legend({'Exponential Fit', 'Simulation', 'Points used for fit'})
xlabel('$t-t_c$ [s]');
ylabel('$a$ [mm]');
xlim(tLimits);
ylim([1e-3, 3e-1]);
fitWindowTicks(tm(fitEnds), a(fitEnds)*1e3, diff(tLimits));

% --- Fig 3b: strain rate --------------------------------------------------
figure();
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
end

function fitWindowTicks(t, a, tSpan)
% short dashed diagonal ticks through the first and last fitted points
dt = 0.02*tSpan;
for k = 1:numel(t)
    plot(t(k) + [-dt, dt], a(k)*[10^-0.25, 10^0.25], 'b--', 'LineWidth', 0.5, ...
        'HandleVisibility', 'off');
end
end
