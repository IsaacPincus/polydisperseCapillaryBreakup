function ax = plotMWDistSplit(x, y, buckets, colours)
%PLOTMWDISTSPLIT Plot a MW distribution with its modes from splitMWDist.
%
%   Shades each bucket and marks its MW with a vertical line whose height
%   is proportional to the bucket's number fraction. colours has one row
%   per bucket (bucket ii uses row N-ii+1).
fig = paperFigure();
ax = gca;
N = numel(buckets.M);
plot(x, y, 'k');
for ii = 1:N
    yii = zeros(size(y));
    range = buckets.edges(ii):buckets.edges(ii+1);
    yii(range) = y(range);
    area(x, yii, 'FaceColor', colours(N-ii+1,:), 'FaceAlpha', 0.2)
end
for ii = 1:N
    plot([buckets.M(ii), buckets.M(ii)], [1e-13, buckets.n(ii)*max(y)/max(buckets.n)], ...
        '-', 'LineWidth', 2, 'Color', colours(N-ii+1,:))
end
xlabel('$M$ [g/mol]')
ylabel('$W$ [mol/g]')
end
