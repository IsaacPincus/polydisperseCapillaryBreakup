function inset = addInset(ax, position, xLimits, yLimits, showTicks)
%ADDINSET Copy axes ax into a zoomed inset drawn on top of it.
inset = copyobj(ax, ax.Parent);
inset.Position = position;
xlabel(inset, '')
ylabel(inset, '')
xlim(inset, xLimits)
ylim(inset, yLimits)
if ~showTicks
    inset.XTick = [];
    inset.YTick = [];
end
end
