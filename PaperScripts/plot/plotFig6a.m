function plotFig6a(r)
%PLOTFIG6A Blend thinning curves for 0, 200 and 500 wppm PS7.
%   The paper's Fig 6a is the 200 wppm PS7 radius plot (second figure).

layout = {  % ylabel, data scale, sim scale, ylim, legend position, textbox position, ep labels
    '$a$ [mm]', 1/2, 1e3, [7e-3 1.005]/2, [0.647 0.543 0.2511 0.2937], [0.6779 0.8335 0.1777 0.07038], 'latex'
    '$a$ [mm]', 1/2, 1e3, [7e-3 1.005]/2, [0.61,0.49,0.25,0.34],       [0.64,0.83,0.3,0.07],          'ppm'
    '$D$ [mm]', 1,   2e3, [7e-3 1.005],   [0.61,0.49,0.25,0.34],       [0.64,0.83,0.3,0.07],          'ppm'
    };
for k = 1:numel(r.blocks)
    b = r.blocks{k};
    [yLabel, dataScale, simScale, yLimits, legendPos, boxPos, epLabels] = layout{k,:};
    colours = winter(numel(b.cPS16)+1);
    labels = concentrationLabels('PS16', b.cPS16);

    plotThinningCurves(num2cell(b.t_exp, 1), num2cell(b.D_exp*dataScale, 1), b.sim, ...
        labels, colours, simScale);
    ylabel(yLabel);
    l1 = legend;
    l1.Position = legendPos;
    xlim([-0.05 0.65]);
    ylim(yLimits);
    annotation('textbox', boxPos, ...
        'String', {sprintf('$c_\\mathrm{PS7} = %d$ [wppm]', b.cPS7(1))}, ...
        'FitBoxToText', 'off');

    if strcmp(epLabels, 'ppm')
        labels = arrayfun(@(c) sprintf('%d ppm PS16', c), b.cPS16, 'UniformOutput', false);
    end
    plotStrainRates(b.tavg_exp, b.ep_exp, b.sim, labels, colours);
end
end
