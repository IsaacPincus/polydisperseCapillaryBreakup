function plotFig2(r)
%PLOTFIG2 SAOS of PS7 (Fig 2a), PS16 (Fig 2b) and their blend (not in paper).

moduliFigure(r.freq, r.PS7, '$c_\mathrm{PS7} = 1000$ [wppm]');
moduliFigure(r.freq, r.PS16, '$c_\mathrm{PS16} = 200$ [wppm]');

% blend: model vs measured, no error bars available
fig = figure();
fig.Position(3:4) = [784.8, 629.6];
hold on
plot(r.freq, sum(r.mix.Gp,1), '-k', 'DisplayName', "$G'$", 'LineWidth', 3)
plot(r.freq, sum(r.mix.Gpp,1), '--k', 'DisplayName', "$G''$", 'LineWidth', 3)
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$\omega$ (1/s)')
ylabel("$G', G''$ (Pa)")
legend
plot(r.mix.exp.freq, r.mix.exp.Gp, 'or', 'DisplayName', "$G'$ exp", ...
    'MarkerFaceColor', 'r')
plot(r.mix.exp.freq, r.mix.exp.Gpp, 'or', 'DisplayName', "$G''$ exp")
annotation('textbox', [0.34,0.83,0.3,0.07], ...
    'String', {'$c_\mathrm{PS16} = 200$ [wppm]', '$c_\mathrm{PS7} = 1000$ [wppm]'}, ...
    'FitBoxToText', 'off');
findfigs
end

function moduliFigure(freq, s, label)
fig = figure();
fig.Position(3:4) = [784.8, 629.6];
hold on
plot(freq, sum(s.Gp,1), '-k', 'DisplayName', "$G'_p$", 'LineWidth', 3)
plot(freq, sum(s.Gpp,1), '--k', 'DisplayName', "$G''_p$", 'LineWidth', 3)
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$\omega$ (1/s)')
ylabel("$G'_p, G''_p$ (Pa)")
legend
errorbar(s.exp.freq, s.exp.Gp, s.exp.Gp_std, 'or', 'DisplayName', "$G'$ exp", ...
    'MarkerFaceColor', 'r')
errorbar(s.exp.freq, s.exp.Gpp, s.exp.Gpp_std, 'or', ...
    'DisplayName', "$G'' - \eta_s \omega$ exp")
annotation('textbox', [0.34,0.83,0.3,0.07], 'String', {label}, 'FitBoxToText', 'off');
xlim([5e-1, 5e2]);
ylim([1e-5, 1e1]);
end
