function [ax, inset] = plotRelaxationRatios(cr16, sim16, exp16, cr7, sim7, exp7, inset16, inset7)
%PLOTRELAXATIONRATIOS tau_EC/tau_Z against c/c* for PS16 and PS7 (Fig 5).
%
%   With inset16, inset7 (f_{i,min} tau_EC/tau_Z of the simulations) an
%   inset shows these against c/c*. Returns the main and inset axes
%   (inset is empty without one).

paperFigure();
ax = gca;
plot(cr16, sim16, 'rs', 'MarkerFaceColor', 'r', 'DisplayName', 'PS16 sim');
plot(cr16, exp16, 'rs', 'DisplayName', 'PS16 exp');
plot(cr7, sim7, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'PS7 sim');
plot(cr7, exp7, 'bo', 'DisplayName', 'PS7 exp');
legend('NumColumns', 2, 'Location', 'northwest')
xlabel('$c/c^*$')
ylabel('$\tau_\mathrm{EC}/\tau_Z$')
set(ax, 'XScale', 'log')
xlim([1e-3, 4e-1]);
ylim([0, ceil(2*max([sim16(:); exp16(:); sim7(:); exp7(:)]))/2]);

inset = [];
if nargin > 6
    inset = axes('Position', [0.24 0.44 0.27 0.29]);
    hold(inset, 'on')
    plot(inset, cr16, inset16, 'rs', 'MarkerFaceColor', 'r');
    plot(inset, cr7, inset7, 'bo', 'MarkerFaceColor', 'b');
    c = [cr16(:); cr7(:)];
    set(inset, 'XScale', 'log')
    xlim(inset, [min(c)/1.5, max(c)*1.5])
    ylim(inset, [0, max(1.5, ceil(2*max([inset16(:); inset7(:)]))/2)])
    xlabel(inset, '$c/c^*$')
    ylabel(inset, '$f_{i,\mathrm{min}}\tau_\mathrm{EC}/\tau_Z$')
end
end
