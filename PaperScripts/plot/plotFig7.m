function plotFig7(r)
%PLOTFIG7 Stress balance (Figs 7a,c) and effective Wi (Figs 7b,d).

% --- Figs 7a,b: one mode per species --------------------------------------
s = r.twoMode;
colours = parula(6);
colours = colours([1,5],:);
tm = s.sol.t - s.relax.t_c;
tLimits = [tm(1), tm(end)];

ax = stressAxes(r.p, s, tm);
for pp = 1:numel(s.modes.M)
    plot(tm, s.sol.stress(:,pp), '-', 'Color', colours(pp,:));
end
plot(tm, viscousInertialStress(r.p, s.sol), 'k-.', 'linewidth', 1.5);
l1 = legend({'Capillary Stress', 'PS16 Polymer Stress', 'PS7 Polymer Stress', ...
    'Viscous and Intertial Stresses'});
xlim(tLimits)
ylim([8e-2, 1e5])
addInset(ax, [0.465454545454546,0.441181102362205,0.218,0.23], ...
    [tm(1), 0.02], [3, 700], false);
l1.Position = [0.162,0.719,0.388,0.174];
annotation(gcf, 'arrow', [0.464426877470356,0.189723320158103], ...
    [0.606882281916595,0.537419634472297]);

paperFigure();
ax = gca;
names = {'PS16', 'PS7'};
for pp = 1:numel(s.modes.M)
    plot(tm, s.sol.Wi(pp,:), '-', 'Color', colours(pp,:), 'DisplayName', names{pp});
end
xlabel('$t-t_c$ [s]')
ylabel('$W\!i_\mathrm{eff}$')
legend('Location', 'northeast')
set(ax, 'YScale', 'log')
xlim(tLimits)
ylim([5e-2, 1e2])
plot(tm, ones(size(tm))*0.5, 'k--', 'DisplayName', '$W\!i_\mathrm{eff} = 0.5$ (coil-stretch)')
addInset(ax, [0.289555555555555,0.654241469816274,0.218,0.23], ...
    [tm(1), 0.01], [2e-1, 9e1], true);
annotation(gcf, 'arrow', [0.280740740740741,0.184782608695652], ...
    [0.801837270341207,0.742782152230971]);

% --- Figs 7c,d: 30 modes ---------------------------------------------------
% thick dashed line: the sub-species that just touches Wi = 0.5
s = r.multiMode;
Nm = numel(s.modes.M);
colours = parula(Nm);
tm = s.sol.t - s.relax.t_c;
tLimits = [tm(1), tm(end)];

stressAxes(r.p, s, tm);
plot(tm, viscousInertialStress(r.p, s.sol), 'k-.', 'linewidth', 1.5);
ylabel('Stress [Pa]')
h = gobjects(1, Nm);
for pp = 1:Nm
    h(pp) = plot(tm, s.sol.stress(:,pp), '-', 'Color', colours(pp,:));
end
legend({'Capillary Stress', 'Viscous and Intertial Stresses'});
xlim(tLimits)
ylim([5e-5, 1e5])
emphasise(h(s.coilStretchMode));

paperFigure();
ax = gca;
for pp = 1:Nm
    h(pp) = plot(tm, s.sol.Wi(pp,:), '-', 'Color', colours(pp,:));
end
xlabel('$t-t_c$ [s]')
ylabel('$W\!i_{\mathrm{eff}, i}$')
set(ax, 'YScale', 'log')
xlim(tLimits)
ylim([1e-2, 1e2])
plot(tm, ones(size(tm))*0.5, 'k--', 'DisplayName', '$W\!i_\mathrm{eff} = 0.5$ (coil-stretch)')
emphasise(h(s.coilStretchMode));
addInset(ax, [0.342888888888889,0.641118110236221,0.218,0.23], ...
    [tm(1), 0.01], [5e-2, 9e1], true);
annotation(gcf, 'arrow', [0.33498023715415,0.184782608695652], ...
    [0.784776902887139,0.742782152230971]);
end

function ax = stressAxes(p, s, tm)
% new figure with the capillary stress chi/a
paperFigure();
ax = gca;
plot(tm, p.chi./s.sol.a, 'k', 'linewidth', 2);
xlabel('$t-t_c$ [s]')
ylabel('stress [Pa]')
set(ax, 'YScale', 'log')
end

function sigma = viscousInertialStress(p, sol)
sigma = 3*p.eta_s*sol.ep + 1/8*p.rho*sol.ep.^2.*sol.a.^2;
end

function emphasise(line)
line.LineStyle = '--';
line.LineWidth = 5;
uistack(line, 'top');
end
