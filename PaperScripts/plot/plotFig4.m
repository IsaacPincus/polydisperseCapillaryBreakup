function plotFig4(r)
%PLOTFIG4 Thinning of single-species solutions, simulation vs experiment
%   (Fig 4a: PS16, Fig 4b: PS7), plus their strain rates (supplementary).

% --- Fig 4a: PS16 ---------------------------------------------------------
s = r.PS16;
colours = winter(numel(s.cPS16)+1);
plotThinningCurves(num2cell(s.t_exp, 1), num2cell(s.D_exp/2, 1), s.sim, ...
    concentrationLabels('PS16', s.cPS16), colours, 1e3);
ylabel('$a$ [mm]');
l1 = legend;
l1.Position = [0.54485,0.5191,0.3131,0.2941];
xlim([-0.05 0.65]);
ylim([5e-3 1.005]/2);
annotation('textbox', [0.586,0.825,0.251,0.070], ...
    'String', {'$c_\mathrm{PS7} = 0$ [wppm]'}, 'FitBoxToText', 'off');

plotStrainRates(s.tavg_exp, s.ep_exp, s.sim, ...
    arrayfun(@(c) sprintf('%d ppm PS16', c), s.cPS16, 'UniformOutput', false), colours);

% --- Fig 4b: PS7 ----------------------------------------------------------
s = r.PS7;
colours = copper(numel(s.cPS7)+1);
colours = flipud(colours(1:end-1,:));
plotThinningCurves(s.t_exp, cellfun(@(D) D/2, s.D_exp, 'UniformOutput', false), ...
    s.sim, concentrationLabels('PS7', s.cPS7), colours, 1e3);
ylabel('$a$ [mm]');
l1 = legend;
l1.Position = [0.551,0.627,0.32,0.19];
xlim([-0.05 0.65]);
ylim([5e-3 1.005]/2);
annotation('textbox', [0.586,0.824,0.30,0.07], ...
    'String', {'$c_\mathrm{PS16} = 0$ [wppm]'}, 'FitBoxToText', 'off');

plotStrainRates(s.tavg_exp, s.ep_exp, s.sim, concentrationLabels('PS7', s.cPS7), colours);
end
