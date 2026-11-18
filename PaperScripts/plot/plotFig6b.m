function plotFig6b(r)
%PLOTFIG6B tau_EC vs PS16 concentration for 0-1000 wppm PS7 (Fig 6b).
%   tau_EC of the simulations (lines) and of the experiments (stars) are both
%   found as in Fig 3: an exponential fit a ~ exp(-t/(3 tau_EC)) in the EC
%   window. The tau_EC/tau_Z plot against c/c* is supplementary.

colours = flip(copper(4));

% --- Fig 6b ----------------------------------------------------------------
fig = paperFigure();
ax = gca;
for k = 1:numel(r.blocks)
    b = r.blocks{k};
    plot(b.cPS16, b.sim.tau_EC, '-', 'Color', colours(k,:));
end
for k = 1:numel(r.blocks)
    b = r.blocks{k};
    plot(b.cPS16, b.tau_EC_exp, 'p', ...
        'Color', colours(k,:), 'MarkerFaceColor', colours(k,:));
end
set(ax, 'XScale', 'log', 'YScale', 'log')
xlabel('$c_\mathrm{PS16}$ [wppm]')
ylabel('$\tau_\mathrm{EC}$ [s]')
legend({' ', ' ', ' ', ' ', ...
    '$\qquad$ 0', '$\qquad$ 200', '$\qquad$ 500', '$\qquad$ 1000'}, ...
    'Position', [0.5379 0.2289 0.30303 0.19773], ...
    'Orientation', 'vertical', 'NumColumns', 2);
legendHeader(fig);

% --- supplementary: tau_EC/tau_Z vs c/c* ----------------------------------
plotRelaxationRatios(r.cr16, r.tau_EC_ratio_PS16, r.tau_EC_exp_ratio_PS16, ...
    r.cr7, r.tau_EC_ratio_PS7, r.tau_EC_exp_ratio_PS7);
end

function legendHeader(fig)
% column headings above the two-column legend
annotation(fig, 'textbox', [0.6408 0.4269 0.0893 0.06356], 'String', {'Exp'});
annotation(fig, 'textbox', [0.5375 0.4249 0.0864 0.063564], 'String', {'Sim'});
annotation(fig, 'textbox', [0.74604 0.42495 0.25 0.06356], 'String', {'$c_\mathrm{PS7}$ [wppm]'});
end
