function plotFig10(r)
%PLOTFIG10 Effect of the plate retraction speed V0.
%   Fig 10a: tau_EC against V0 (same y-limits as Fig 9), with the legend.
%   Fig 10b: polymer stress of each species when the plates stop, divided
%            by the mid-filament strain rate then, 1.5 V0/L_end (Appendix
%            B), i.e. an effective polymer extensional viscosity, against
%            V0. PS16: solid lines and filled circles, coloured by c_PS16;
%            PS7: black dashed lines and open squares, the only legend entry.

V0 = r.retraction.V0_vals*1000;   % mm/s

% --- Fig 10a --------------------------------------------------------------
tauEcFigure(r, V0, r.retraction.tau_EC, '$V_0$ [mm/s]');
lg = legend('Location', 'southeast', 'NumColumns', 2);
% c_PS7 label just above the legend, right-aligned with it
w = 0.4;
annotation('textbox', [sum(lg.Position([1 3]))-w, sum(lg.Position([2 4]))+0.01, w, 0.07], ...
    'String', sprintf('$c_\\mathrm{PS7} = %g$ [wppm]', r.cPS7), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')

% --- Fig 10b --------------------------------------------------------------
paperFigure();
ax = gca;
colours = winter(numel(r.cPS16));
rate = r.retraction.epend;   % mid-filament strain rate when the plates stop [1/s]
for k = 1:numel(r.cPS16)
    plot(V0, r.retraction.stress16(:,k)./rate(:,k), '-o', ...
        'Color', colours(k,:), 'MarkerFaceColor', colours(k,:), 'HandleVisibility', 'off');
    h7 = plot(V0, r.retraction.stress7(:,k)./rate(:,k), '--sk', 'HandleVisibility', 'off');
end
set(h7, 'HandleVisibility', 'on', ...   % one legend entry for PS7
    'DisplayName', sprintf('PS7, $c_\\mathrm{PS7} = %g$ [wppm]', r.cPS7));
set(ax, 'XScale', 'log', 'YScale', 'log', 'XMinorTick', 'on')
xlabel('$V_0$ [mm/s]')
ylabel('$\Delta \sigma_{p,i}/(V_0/L_\mathrm{end})$ at end of pre-stretch [Pa s]')
legend('Location', 'best')
end
