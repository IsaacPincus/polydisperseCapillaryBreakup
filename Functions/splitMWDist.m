function [n, M, mol_polymer_density] = ...
    splitMWDist(graphs_on, colours,...
    x, y, noPoints, maxMW, N, nu, density_mass_PS)

    % given a molecular weight distribution (weight-fraction vs MW), splits
    % the distribution into N buckets, and finds the viscosity-averaged MW
    % of each bucket, as well as the number density of polymers within that
    % bucket (normalised to sum to 1).
    % Inputs:
    %   graphs_on - plots MW distribution if true
    %   colours - colours to use for graphs
    %   x, y - weight-fraction probability density (y), MW (x)
    %   noPoints - Number of points in MW dist [length(x)]
    %   maxMW - not necessary for this function
    %   N - number of buckets to split MW dist into
    %   nu - exponent in viscosity-averaged MW
    %   density_mass_PS - total mass density of polymer in solution
    
    mol_polymer_density = trapz(x, y*density_mass_PS./x);

    % instead just subdivide regions
    if graphs_on
        figure10 = figure();
        figure10.Position = [488,152.2,810,610];
        axes1 = axes('Parent',figure10);
        hold(axes1,'on');
        plot(x,y, 'k');
    end
    
    indices = round(linspace(1,noPoints, N+1));
    
    xvals = x(indices);
    yvals = y(indices);
    
    for ii = 1:N
        yii = zeros(size(y));
        N_range = (indices(ii):indices(ii+1));
        yii(N_range) = y(N_range);
%         % viscosity-averaged MW
%         ML(ii) = (trapz(x, yii.*x.^(nu))/trapz(x, yii))^(1/nu);
        % extensibility-averaged MW
        ML(ii) = (trapz(x, yii.*x.^(1+nu))/trapz(x, yii))^(1/(1+nu));
        if isnan(ML(ii))
            ML(ii) = mean(x(N_range));
        end
        ni(ii) = trapz(x, yii./x);
        colour = colours(N-ii+1,:);
        if graphs_on
            area(x,yii, 'FaceColor',colour, 'FaceAlpha',0.2)
        end
    end
    total_n = trapz(x,y./x);
    n = ni/total_n;
    M = ML;

    for ii = 1:N
        if graphs_on
            colour = colours(N-ii+1,:);
            plot([ML(ii),ML(ii)], [1e-13,n(ii)*max(y)/max(n)], ...
                '-', 'LineWidth',2,'Color',colour)
        end
    end

    % sort modes so they are in descending order by MW
    [M, indices] = sort(M, 'descend');
    n = n(indices);
    
    if graphs_on
        xlabel('$M$ [g/mol]')
        ylabel('$W$ [mol/g]')
    %     axes1.YScale = 'log';
        
        hold(axes1,'off');
    end
   

end


