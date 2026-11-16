function plotFig6b(r)
%PLOTFIG6B tau_EC vs PS16 concentration for 0-1000 wppm PS7.
%   The paper's Fig 6b is the first figure (tau_EC from min(ep), with the
%   Haward et al. data). The others are supplementary.

colours = flip(copper(4));

% --- Fig 6b: tau_EC from the strain-rate minimum --------------------------
fig58 = paperFigure();
ax58 = gca;
for k = 1:numel(r.blocks)
    b = r.blocks{k};
    plot(b.cPS16, b.sim.tau_ep_min, '-', 'Color', colours(k,:));
end
for k = 1:4
    plot(r.haward(:,k*2-1), r.haward(:,k*2)/1000, 'p', ...
        'Color', colours(k,:), 'MarkerFaceColor', colours(k,:));
end
set(ax58, 'XScale', 'log', 'YScale', 'log')
xlabel('$c_\mathrm{PS16}$ [wppm]')
ylabel('$\tau_\mathrm{EC}$ [s]')
legend({' ', ' ', ' ', ' ', ...
    '$\qquad$ 0', '$\qquad$ 200', '$\qquad$ 500', '$\qquad$ 1000'}, ...
    'Position', [0.5379 0.2289 0.30303 0.19773], ...
    'Orientation', 'vertical', 'NumColumns', 2);
legendHeader(fig58, 0, '$\tau_\mathrm{EC}$ from $\min(\dot{\epsilon})$');

% --- supplementary: tau_EC from exponential fits --------------------------
fig59 = paperFigure();
ax59 = gca;
for k = 1:numel(r.blocks)
    b = r.blocks{k};
    plot(b.cPS16, b.sim.tau_EC, ':', 'Color', colours(k,:));
    plot(b.cPS16, b.tau_EC_exp, 'o', 'Color', colours(k,:));
end
set(ax59, 'XScale', 'log', 'YScale', 'log')
xlabel('$c_\mathrm{PS16}$ [wppm]')
ylabel('$\tau_\mathrm{EC}$ [s]')
legend({' ', '$\qquad$ 0', ' ', '$\qquad$ 200', ' ', '$\qquad$ 500', ...
    ' ', '$\qquad$ 1000'}, ...
    'Position', [0.5579 0.2289 0.30303 0.19773], ...
    'Orientation', 'horizontal', 'NumColumns', 2);
legendHeader(fig59, 0.02, '$\tau_\mathrm{EC}$ from $D \sim e^{-t/3\tau_\mathrm{EC}}$');

% --- supplementary: both together -----------------------------------------
figure();
hold on
copyobj(ax58.Children, gca);
copyobj(ax59.Children, gca);
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$c_\mathrm{PS16}$ [wppm]')
ylabel('$\tau_\mathrm{EC}$ [s]')

% --- supplementary: tau_EC/tau_Z vs c/c* ----------------------------------
plotRelaxationRatios(r.cr16, r.tau_EC_ratio_PS16, r.tau_EC_exp_ratio_PS16, ...
    r.cr7, r.tau_EC_ratio_PS7, r.tau_EC_exp_ratio_PS7);
end

function legendHeader(fig, dx, caption)
% column headings above the two-column legend
annotation(fig, 'textbox', [0.6408+dx 0.4269 0.0893 0.06356], 'String', {'Exp'});
annotation(fig, 'textbox', [0.5375+dx 0.4249 0.0864 0.063564], 'String', {'Sim'});
annotation(fig, 'textbox', [0.74604+dx 0.42495 0.25 0.06356], 'String', {'$c_\mathrm{PS7}$ [wppm]'});
annotation(fig, 'textbox', [0.54604+dx 0.15495 0.5 0.06356], 'String', {caption});
end
