function plotFig1(r)
%PLOTFIG1 Fig 1 (PS7 and PS16 as one mode each) and Fig A1 (30 modes).

% --- Fig 1 ----------------------------------------------------------------
colours = parula(6);
cPS16 = colours(1,:);
cPS7 = colours(5,:);
m = r.modes;
heightPS7 = max(r.yPS7);
heightPS16 = heightPS7*m.n(1)/m.n(2);   % number ratio PS16:PS7 in the blend

paperFigure();
plot(r.x, r.yPS7, 'k');
area(r.x, r.yPS7, 'FaceColor', cPS7, 'FaceAlpha', 0.2)
plot([m.MPS7, m.MPS7], [1e-13, heightPS7], '-', 'LineWidth', 2, 'Color', cPS7)
plot(r.x, r.yPS16, 'k');
area(r.x, r.yPS16, 'FaceColor', cPS16, 'FaceAlpha', 0.2)
plot([m.MPS16, m.MPS16], [0, heightPS16], '-', 'LineWidth', 2, 'Color', cPS16)
xlabel('$M$ [g/mol]')
ylabel('$W$ [mol/g]')

annotation(gcf, 'textbox', [0.361660079051385,0.782152230971129, ...
    0.0800395256917,0.05748031496063], 'String', {'PS7'});
annotation(gcf, 'textarrow', [0.7549,0.59930], [0.28740,0.21378], ...
    'String', {'$M_{L,i=1}$'})
annotation(gcf, 'textbox', [0.65316205533597,0.568241469816274, ...
    0.0899209486166,0.058792650918636], 'String', {'PS16'})
annotation(gcf, 'textarrow', [0.4457,0.3454], [0.60412,0.54055], ...
    'String', {'$M_{L,i=2}$'})
annotation(gcf, 'textbox', [0.586956,0.769028,0.27470355,0.11653543], ...
    'String', {sprintf('$c_\\mathrm{PS7} = %d$ [wppm]', r.cPS7), ...
    sprintf('$c_\\mathrm{PS16} = %d$ [wppm]', r.cPS16)});

% --- Fig A1 ---------------------------------------------------------------
J = numel(r.A1.M);
plotMWDistSplit(r.x, r.A1.y, r.A1.buckets, parula(J));
annotation(gcf, 'textbox', [0.5 0.75 0.4 0.15], ...
    'String', {sprintf('$c_\\mathrm{PS7} = %d$ [wppm]', r.cPS7), ...
    sprintf('$c_\\mathrm{PS16} = %d$ [wppm]', r.cPS16), sprintf('$J = %d$', J)});
annotation(gcf, 'textarrow', [0.472332015810277 0.350790513833992], ...
    [0.583224115334207 0.516382699868938], 'String', {'number density'});
annotation(gcf, 'textarrow', [0.563241106719368 0.596837944664032], ...
    [0.397116644823067 0.292267365661861], 'String', {'mass density'})
end
