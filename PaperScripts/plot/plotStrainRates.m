function plotStrainRates(tavg_exp, ep_exp, sim, labels, colours)
%PLOTSTRAINRATES Measured (symbols) and simulated (lines) strain rates.
paperFigure();
for ii = 1:numel(ep_exp)
    plot(tavg_exp{ii}, ep_exp{ii}, 's', 'Color', colours(ii,:), 'DisplayName', labels{ii})
    plot(sim.tModel{ii} - sim.t_c(ii), sim.epModel{ii}, '-', 'linewidth', 2, ...
        'Color', colours(ii,:), 'HandleVisibility', 'off');
end
set(gca, 'YScale', 'log')
xlabel('$t-t_c$ [s]');
ylabel('$\dot{\epsilon}$ [1/s]');
xlim([-0.05 0.65]);
ylim([4, 2000]);
end
