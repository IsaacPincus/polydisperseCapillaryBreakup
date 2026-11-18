function plotFigB1(r)
%PLOTFIGB1 Liquid-bridge shape during the plate pre-stretch (Appendix B).
%   (a) Free surface of the bridge, initially a cylinder of radius R0
%       between plates L0 apart, at a few plate gaps up to the end of the
%       pre-stretch (mid-filament radius a0). Circles: Lagrangian points on
%       the surface, open at the start and filled at the end, with their
%       paths (dotted); R0 and L0 mark the initial bridge. The plates are
%       drawn as rectangles of radius R0, so the contact lines stay
%       pinned at their edges; earlier positions of the top plate are
%       outlined in the grey of the matching surface, labelled with the
%       time at which it is reached for V0 = r.V0.
%   (b) z0 (Eq. B7), 1/(1-4z0/H), 4z0/(H-4z0) and (1-4z0/H)/(1-4z0) against
%       the starting height zeta0, at the end of the pre-stretch.
%   (c) The same information as (a) drawn as a diagram: only the z axis,
%       the final bridge filled, the initial cylinder dashed and the axis of
%       symmetry dash-dotted.

% --- (a) bridge shape -----------------------------------------------------
fig = paperFigure();
ax = gca;
mm = 1e3;
nP = numel(r.H_profiles);
shades = flip(gray(nP + 2));
shades = shades(3:end, :);   % light grey to black, ending at the final shape
R0 = r.R0*mm;
t = 0.3;   % plate thickness drawn [mm]
for k = 1:nP
    rk = r.profile.r(:,k)*mm;
    zk = r.profile.z(:,k)*mm;
    plot([rk; NaN; -rk], [zk; NaN; zk], '-', 'Color', shades(k,:), ...
        'HandleVisibility', 'off');
    if k < nP   % earlier top-plate positions: outline only
        plate(R0, max(zk), t, 'none', shades(k,:));
    end
    text(R0 + 0.2, max(zk) + t/2, sprintf('$t = %.3g$ s', r.t_profiles(k)), ...
        'Color', shades(k,:), 'VerticalAlignment', 'middle');
end
plate(R0, 0, -t, 0.7*[1 1 1], 'k');                   % bottom plate
plate(R0, r.Lend*mm, t, 0.7*[1 1 1], 'k');            % top plate at the end
plot([r.path.r; NaN(1, numel(r.zeta0))]*mm, [r.path.z; NaN(1, numel(r.zeta0))]*mm, ...
    'k:', 'HandleVisibility', 'off');
plot(r.path.r(1,:)*mm, r.path.z(1,:)*mm, 'ko', 'DisplayName', 'Initial ($H = 1$)');
plot(r.path.r(end,:)*mm, r.path.z(end,:)*mm, 'ko', 'MarkerFaceColor', 'k', ...
    'DisplayName', sprintf('Final ($H = %.3g$)', r.Lend/r.L0));
% mid-filament radius a0 at the end
zMid = r.Lend/2*mm;
plot(r.a0*mm*[-1, 1], zMid*[1, 1], 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
text(0, zMid, '$2a_0$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
axis equal
% directions of the dimensionless coordinates xi = r/R0 and zeta = z/L0
% and of the velocities u and v, from the axis on the bottom plate
arrowLength = 0.8;   % mm
quiver(0, 0, arrowLength, 0, 0, 'k', 'MaxHeadSize', 0.5, 'HandleVisibility', 'off');
quiver(0, 0, 0, arrowLength, 0, 'k', 'MaxHeadSize', 0.5, 'HandleVisibility', 'off');
text(arrowLength + 0.05, 0, '$\xi, u$', 'VerticalAlignment', 'bottom');
text(0.05, arrowLength + 0.05, '$\zeta, v$', 'HorizontalAlignment', 'left');
% dimensions of the initial bridge: R0 below the bottom plate, L0 at the left
L0 = r.L0*mm;
dimension([0, R0], -t - 0.35*[1, 1], '$R_0$', 'center', 'top');
dimension(-R0 - 0.35*[1, 1], [0, L0], '$L_0$', 'right', 'middle');
xlim(ax, R0*[-1.6, 1.9])   % room for L0 on the left and the time labels on the right
ylim(ax, [-t - 1, r.Lend*mm + t + 0.2])
% narrow the figure so the equal-axis plot fills it (keeps its height)
fig.Position(3) = fig.Position(4)*ax.Position(4)/ax.Position(3) ...
    *diff(ax.XLim)/diff(ax.YLim);
xlabel('$r$ [mm]')
ylabel('$z$ [mm]')
set(findobj(ax, 'Type', 'text'), 'FontSize', 16)   % labels inside the plot

plotZ0(r);
bridgeDiagram(r);
end

function bridgeDiagram(r)
% (c) diagram version of (a)
fig = paperFigure();
ax = gca;
mm = 1e3;
R0 = r.R0*mm;
L0 = r.L0*mm;
Lend = r.Lend*mm;
t = 0.3;   % plate thickness drawn [mm]
nP = numel(r.H_profiles);
shades = flip(gray(nP + 2));
shades = shades(3:end, :);

% final bridge, filled
rEnd = r.profile.r(:,end)*mm;
zEnd = r.profile.z(:,end)*mm;
fillColour = [0.78 0.88 1];
patch([rEnd; -flipud(rEnd)], [zEnd; flipud(zEnd)], fillColour, ...
    'EdgeColor', 'k', 'HandleVisibility', 'off');
% initial cylinder, dashed
% plot([-R0, -R0, R0, R0], [0, L0, L0, 0], 'k--', 'LineWidth', 0.75, ...
%     'HandleVisibility', 'off');
plate(R0, L0, t, 'none', shades(1,:));   % top plate at t = 0
% intermediate surfaces and top-plate positions, with their times
for k = 1:nP-1
    rk = r.profile.r(:,k)*mm;
    zk = r.profile.z(:,k)*mm;
    plot([rk; NaN; -rk], [zk; NaN; zk], '-', 'Color', shades(k,:), ...
        'HandleVisibility', 'off');
    plate(R0, max(zk), t, 'none', shades(k,:));
end
for k = 1:nP
    if k == nP
        text(R0 + 0.2, max(r.profile.z(:,k))*mm + t/2, ...
            sprintf('$t = t_\\mathrm{end}$ \n$ = %.3g$ s',...
            r.t_profiles(k)), 'Color', shades(k,:), ...
            'VerticalAlignment', 'middle');
    else
        text(R0 + 0.2, max(r.profile.z(:,k))*mm + t/2, ...
            sprintf('$t = %.3g$ s', r.t_profiles(k)), 'Color', shades(k,:), ...
            'VerticalAlignment', 'middle');
    end
end
plate(R0, 0, -t, 0.7*[1 1 1], 'k');
plate(R0, Lend, t, 0.7*[1 1 1], 'k');
% axis of symmetry
plot([0, 0], [-t - 0.2, Lend + t + 0.2], 'k-.', 'LineWidth', 0.5, ...
    'HandleVisibility', 'off');
% Lagrangian points and their paths; the midplane point (zeta0 = 0.5) as
% diamonds on a darker dotted path, the others fainter
mid = abs(r.zeta0 - 0.5) < 1e-9;
other = ~mid;
plot([r.path.r(:,other); NaN(1, nnz(other))]*mm, ...
    [r.path.z(:,other); NaN(1, nnz(other))]*mm, ...
    ':', 'Color', 0.45*[1 1 1], 'LineWidth', 0.75, 'HandleVisibility', 'off');
plot(r.path.r(1,other)*mm, r.path.z(1,other)*mm, 'ko', 'MarkerFaceColor', 'w', ...
    'HandleVisibility', 'off');
plot(r.path.r(end,other)*mm, r.path.z(end,other)*mm, 'ko', 'MarkerFaceColor', 'k', ...
    'HandleVisibility', 'off');
plot(r.path.r(:,mid)*mm, r.path.z(:,mid)*mm, 'k:', 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');
plot(r.path.r(1,mid)*mm, r.path.z(1,mid)*mm, 'kd', 'MarkerSize', 9, ...
    'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
plot(r.path.r(end,mid)*mm, r.path.z(end,mid)*mm, 'kd', 'MarkerSize', 9, ...
    'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
% 2a0 at the midplane, R0 and L0 of the initial bridge
zMid = Lend/2;
% double-headed arrow across the neck, from the centreline out to each side
quiver([0, 0], zMid*[1, 1], r.a0*mm*[-1, 1], [0, 0], 0, 'k', 'LineWidth', 0.5, ...
    'MaxHeadSize', 0.4, 'HandleVisibility', 'off');
text(0, zMid, '$2a_0$', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'BackgroundColor', fillColour, 'Margin', 1);
dimension([0, R0], -t - 0.35*[1, 1], '$R_0$', 'center', 'top');
dimension(-R0 - 0.35*[1, 1], [0, L0], '$L_0$', 'right', 'middle');
dimension(-R0 - 1.1*[1, 1], [0, Lend], '$L_\mathrm{end}$', 'right', 'middle');
% coordinate directions from the axis on the bottom plate
arrowLength = 0.8;
quiver(0, 0, arrowLength, 0, 0, 'k', 'MaxHeadSize', 0.5, 'HandleVisibility', 'off');
quiver(0, 0, 0, arrowLength, 0, 'k', 'MaxHeadSize', 0.5, 'HandleVisibility', 'off');
text(arrowLength + 0.05, 0, '$\xi, u$', 'VerticalAlignment', 'bottom');
text(0.05, arrowLength + 0.05, '$\zeta, v$', 'HorizontalAlignment', 'left');

axis equal
box(ax, 'off')
% ax.XAxis.Visible = 'off';   % only the z axis, on the left
xlabel('$r$ [mm]')
ylabel('$z$ [mm]')
xlim(ax, R0*[-1.85, 1.9])   % room for L0 and Lend on the left
ylim(ax, [-t - 0.9, Lend + t + 0.3])
fig.Position(3) = fig.Position(4)*ax.Position(4)/ax.Position(3) ...
    *diff(ax.XLim)/diff(ax.YLim);
set(findobj(ax, 'Type', 'text'), 'FontSize', 16)
end

function dimension(x, z, label, hAlign, vAlign)
% dimension line from (x(1), z(1)) to (x(2), z(2)) with end ticks, labelled
% beside its middle (below it for a horizontal line, left for a vertical one)
tick = 0.12;
plot(x, z, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
if z(1) == z(2)   % horizontal
    plot([x; x], z(1) + tick*[-1; 1]*[1, 1], 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    text(mean(x), z(1) - tick, label, 'HorizontalAlignment', hAlign, 'VerticalAlignment', vAlign);
else              % vertical
    plot(x(1) + tick*[-1; 1]*[1, 1], [z; z], 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    text(x(1) - tick, mean(z), label, 'HorizontalAlignment', hAlign, 'VerticalAlignment', vAlign);
end
end

function plate(R0, zFace, thickness, faceColour, edgeColour)
% plate of radius R0 whose wetted face is at z = zFace; it extends a
% distance thickness away from the liquid (negative: below the face)
z = zFace + [0, thickness];
patch([-R0, R0, R0, -R0], z([1 1 2 2]), 'k', 'FaceColor', faceColour, ...
    'EdgeColor', edgeColour, 'HandleVisibility', 'off');
end

function plotZ0(r)
% --- (b) z0 and the combinations of it in Eqs. B6-B10, against zeta0 -----
% at the end of the pre-stretch (H = Lend/L0)
paperFigure();
ax = gca;
H = r.Lend/r.L0;
% H = 1;
zeta0 = r.z0curve.zeta0;
z0 = r.z0curve.z0;
plot(-z0, zeta0, 'k-', 'DisplayName', '$-z_0$');
plot(1./(1 - 4*z0/H), zeta0,  'k--', ...
    'DisplayName', '$1/(1-4z_0/H)$');
% plot(zeta0, 4*z0./(H - 4*z0), 'k:', 'DisplayName', '$4z_0/(H-4z_0)$');
plot((1 - 4*z0/H)./(1 - 4*z0), zeta0,  'k-.', ...
    'DisplayName', '$(1-4z_0/H)/(1-4z_0)$');
ylim(ax, [0, 1])
ylabel('$\zeta_0$')
legend('Location', 'best', 'NumColumns', 1)
text(0.03, 1, sprintf('$H = L_\\mathrm{end}/L_0 = %.3g$', H), ...
    'VerticalAlignment', 'top', 'FontSize', 16);
% set(gca, 'XScale', 'log')
% xlim([1e-3, 1e1])
xlim([0, 2])
end
