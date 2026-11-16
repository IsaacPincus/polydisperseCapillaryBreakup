function plotFig5(r)
%PLOTFIG5 tau_EC/tau_Z against c/c* for PS16 and PS7 (Fig 5), with dashed
%   lines at 1 and an inset of f_{i,min} tau_EC/tau_Z (f_{i,min}: Eq. 22 for
%   species i at the strain-rate minimum).

[ax, inset] = plotRelaxationRatios(r.cr16, r.PS16.tau_EC_ratio, r.PS16.tau_EC_exp_ratio, ...
    r.cr7, r.PS7.tau_EC_ratio, r.PS7.tau_EC_exp_ratio, ...
    r.PS16.f_min_ratio, r.PS7.f_min_ratio);
set(findobj(ax, 'Type', 'line'), 'MarkerSize', 10);
yline(ax, 1, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
yline(inset, 1, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
end
