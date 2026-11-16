function plotFig8(r)
%PLOTFIG8 Weighted mean MW vs c_PS16 (Fig 8b), the corresponding Zimm
%   times (supplementary) and weighted MW distributions (Fig 8a).
%   'terminal' is the weighting by Eq. 33, sigma_p = 2 g tau L^2 ep.

p = r.p;
w = r.sweep;

% --- Fig 8b: mean MW ------------------------------------------------------
paperFigure();
weightedLines(r.cPS16, w.meanStressP, w.meanStressAna, w.meanNumber);
xlabel('$c_\mathrm{PS16}$ [ppm]')
ylabel('$M$ [MDa]')
allM = [w.meanStressP(:); w.meanStressAna(:); w.meanNumber(:)];
ylim([floor(min(allM)), ceil(max(allM))])
sweepLegend(r.cPS7);

% --- supplementary: Zimm time of the mean MW vs tau_EC ----------------------
zimmTime = @(M_MDa) 1/p.S*p.K*(M_MDa*1e6).^(1 + p.alpha)*p.eta_s/(p.NA*p.kB*p.T);   % Eq. 17
paperFigure();
weightedLines(r.cPS16, zimmTime(w.meanStressP), zimmTime(w.meanStressAna), ...
    zimmTime(w.meanNumber));
plot(r.cPS16, w.tau_EC, 'x', 'Color', [0.5 0.5 0.5], 'DisplayName', '$\tau_\mathrm{EC}$')
xlabel('$c_\mathrm{PS16}$ [ppm]')
ylabel('$\tau_Z(M)$ [s]')
sweepLegend(r.cPS7);

% --- Fig 8a ---------------------------------------------------------------
plotWeightedDistributions(r.single);
end

function weightedLines(c, stress, terminal, number)
% one colour per PS7 concentration: stress-weighted with markers,
% terminal-weighted dotted, number-weighted dashed
colours = flip(copper(size(stress, 2)));
for k = 1:size(stress, 2)
    plot(c, stress(:,k), '-o', 'Color', colours(k,:), 'MarkerSize', 4, ...
        'MarkerFaceColor', colours(k,:))
end
for k = 1:size(stress, 2)
    plot(c, terminal(:,k), ':', 'Color', colours(k,:))
end
for k = 1:size(stress, 2)
    plot(c, number(:,k), '--', 'Color', colours(k,:))
end
end

function sweepLegend(cPS7)
% three legend columns (sigma_p, terminal, n), labelled by c_PS7 on the right
n = numel(cPS7);
entries = [repmat({' '}, 1, 2*n), arrayfun(@(c) sprintf('$\\qquad$%d', c), cPS7, ...
    'UniformOutput', false)];
legend(entries, 'Position', [0.4327 0.29017 0.40472 0.1499], 'NumColumns', 3)
annotation('textbox', [0.42888 0.43309 0.46 0.066229], 'String', ...
    ['$\sigma_p$ $\quad$ $\sigma_p^\mathrm{terminal}$ $\quad$ $n$' ...
    ' $\quad$ $c_\mathrm{PS7}$ [ppm]'], ...
    'FitBoxToText', 'off', 'BackgroundColor', [1,1,1]);
end

function plotWeightedDistributions(s)
fig = figure;
fig.Position = [488, 152.2, 810, 610];

ax1 = axes('Parent', fig);
hold(ax1, 'on');
h = plot(s.X1, s.YMatrix1, 'Parent', ax1);
set(h(1), 'DisplayName', '$\sigma_p^\mathrm{terminal}$-weighted', 'Color', [0 0 0], 'LineStyle', ':');
set(h(2), 'DisplayName', '$\sigma_p$-weighted', 'Color', [0 0 0], 'LineStyle', '-', ...
    'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', [0 0 0]);
set(h(3), 'DisplayName', '$n$-weighted', 'Color', [0 0 0], 'LineStyle', '--');
ylim(ax1, [0, 1.05*max(s.YMatrix1(:))]);
ylabel('pdf');
xlabel('$M$ [MDa]');
hold(ax1, 'off');
legend1 = legend(ax1, 'show');
set(legend1, 'Position', [0.145169485784945,0.756049340785143,0.253417204319446,0.149541631316831]);

% inset: capillary, viscous + inertial and modal polymer stresses, with the
% time of the strain-rate minimum; limits follow the simulation
ax2 = axes('Parent', fig, ...
    'Position', [0.531715463816916,0.597808865966411,0.281561264822132,0.298693750105661]);
hold(ax2, 'on');
ax2.FontSize = 13;
h = semilogy(s.X2, s.YMatrix2);
set(h(1), 'DisplayName', 'Capillary', 'LineWidth', 2, 'Color', [0 0 0]);
set(h(2), 'DisplayName', 'Viscous \& Inertial', 'LineStyle', '--', 'Color', [0 0 0]);
modeColours = parula(numel(h) - 2);
for k = 3:numel(h)
    set(h(k), 'Color', modeColours(k-2,:));
end
% y range: every mode's stress at the strain-rate minimum up to the capillary stress
k = find(s.X2 >= s.X3(1), 1);
modal = s.YMatrix2(k, 3:end);
yLimits = [10^floor(log10(min(modal(modal > 0)))), 10^ceil(log10(max(s.YMatrix2(:,1))))];
semilogy(s.X3, yLimits, 'LineWidth', 2, 'Color', [1 0 0]);
text(ax2, s.X3(1), sqrt(yLimits(1)*yLimits(2)), ' $t(\dot{\epsilon}_\mathrm{min})$', ...
    'Color', [1 0 0], 'FontSize', 13);
ylabel('Stress [Pa]');
xlabel('$t-t_c$ [s]');
xlim(ax2, [s.X2(1), s.X2(end)]);
ylim(ax2, yLimits);
hold(ax2, 'off');
set(ax2, 'YMinorTick', 'on', 'YScale', 'log');
end
