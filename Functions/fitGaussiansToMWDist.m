function [n, M, mol_polymer_density] = ...
    fitGaussiansToMWDist(graphs_on, colours,...
    x, y, noPoints, maxMW, N, nu, density_mass_PS)

    if graphs_on
        figure(40);
        hold on
        xlim([0,1]);
    end
    
    disp('hello')
    [x,y] = prepareCurveData(x, y);
    disp('hello')
    x = x/maxMW;
    y = y*maxMW;
    % plot(x,y)
    
    [coeffs, str] = multipleGaussianString(N);
    disp(str)
    disp(coeffs)
    
    myfittype = fittype(str, ...
        dependent="y", independent="x", coefficients = coeffs);
    initial_vals = get_start_points(N,x,y);
    opts = fitoptions('StartPoint', initial_vals, ...
        'Method', 'NonlinearLeastSquares',...
        'TolFun', 1e-12, 'TolX', 1e-12,...
        'Lower', zeros(3*N-2, 1),...
        'Weights', []);
    myfit = fit(x, y, myfittype, opts);
    if graphs_on
        plot(myfit,x,y)
    end
    [a,b,c] = getVars(N,myfit);
    denom = sum(a.*c)*sqrt(2*pi);
    for ii = 1:N
        yii = a(ii)*exp(-(x-b(ii)).^2/(2*c(ii)^2))/denom;
        ML(ii) = (trapz(x, yii.*x.^(nu))/trapz(x, yii))^(1/nu);
        ni(ii) = trapz(x, yii./x);
    end
    total_n = trapz(x,y./x);
    n = ni/total_n;
    M = ML*maxMW;
    mol_polymer_density = trapz(x*maxMW, y/maxMW*density_mass_PS./(x*maxMW));
    
    % sort modes so they are in ascending order by MW
    [M, indices] = sort(M, 'descend');
    n = n(indices);
    % colours = colours(indices,:);
    
    for ii=1:N
        yii = a(ii)*exp(-(x-b(ii)).^2/(2*c(ii)^2))/denom;
        if graphs_on
            plot(x,yii, ':', 'color', colours(N-ii+1,:), ...
                'HandleVisibility', 'off');
        end
    end
    
    if graphs_on
        xlabel('Normalized MW')
        ylabel('Weight-fraction PDF')
    end
    
    if graphs_on
        figure(30);
        hold on
        xlim([0,1]);
    end
    sum_yii = zeros(noPoints,1);
    for ii=1:N
        yii = a(ii)*exp(-(x-b(ii)).^2/(2*c(ii)^2))/denom;
        if graphs_on
            plot(x,yii./x, ':', 'color', colours(N-ii+1,:), ...
                'HandleVisibility', 'off');
        end
        sum_yii = sum_yii + yii./x;
    end
    if graphs_on
        plot(x, sum_yii, 'r-', 'LineWidth',2);
        
        xlabel('Normalized MW')
        ylabel('Number-fraction PDF')
    end
end


%% functions for dist fitting
function [coeffs, str] = multipleGaussianString(N)
% creates a string representing a sum of Gaussians with normalised area

for ii = 1:N
    if ii == 1
        str_upper = "exp(-(x-b1)^2/(2*c1^2))";
        str_lower = "c1";
        coeffs = ["b1", "c1"];
    else
        str_upper = str_upper ...
            + sprintf(" + a%d*exp(-(x-b%d)^2/(2*c%d^2))", ii,ii,ii);
        str_lower = str_lower ...
            + sprintf(" + a%d*c%d", ii, ii);
        coeffs = [coeffs, ...
            sprintf("a%d", ii), ...
            sprintf("b%d", ii), ...
            sprintf("c%d", ii)];
    end
end

str = "(" + str_upper + ")/((" + str_lower + ")*sqrt(2*pi))";

end

function initial = get_start_points(N,x,y)
% choose start points as guess for fitting procedure
lenx = length(x);

% split the interval into N+1 parts, so we have N+2 'cuts' including the
% end points
indices = round(linspace(1,lenx, N+2));

xvals = x(indices);
yvals = y(indices);

bvals = xvals(2:N+1);
avals = yvals(2:N+1);
cval = mean(diff(xvals))/3;

initial = [bvals(1), cval];
count = 3;
for ii = 2:N
    initial(count:count+2) = [avals(ii)/avals(1), bvals(ii), cval];
    count = count + 3;
end

end

% function denom = getDenom(N, fit)
%     c(1) = fit.c1;
%     a(1) = 1;
%     for ii = 2:N
%         c(ii) = fit.(sprintf("c%d", ii));
%         a(ii) = fit.(sprintf("a%d", ii));
%     end
%     denom = sum(a.*c);
% end

function [a,b,c] = getVars(N,fit)
    a(1) = 1;    
    b(1) = fit.b1;
    c(1) = fit.c1;

    for ii = 2:N
        a(ii) = fit.(sprintf("a%d", ii));
        b(ii) = fit.(sprintf("b%d", ii));
        c(ii) = fit.(sprintf("c%d", ii));
    end
end
