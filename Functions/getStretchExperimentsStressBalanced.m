function [t0, y0, a0] = getStretchExperimentsStressBalanced(...
        chi, X, mu, rho,Nm, Ns, L, lambda, g, texp, Dexp)
    % getStretchExperiments runs the FENE-P model with a fixed rate of
    % extension set by experimental data. Once the experimental data shows
    % a maximum in the extension rate, the filament radius and state of
    % extension is output for the 0D filament thinning model.

    % get extension rates
    texp = texp(~isnan(texp));
    Dexp = Dexp(~isnan(texp));
%     figure();
%     hold on
%     plot(texp,Dexp, 'o');
    % smooth out data with spline interpolation
%     [~,DexpSmooth] = spaps(texp,Dexp,3e-13);
%     DexpSmooth = DexpSmooth';
%     [pp,P] = csaps(texp,Dexp);
%     DexpSmooth = fnval(csaps(texp,Dexp,P),texp);
%     plot(texp,DexpSmooth)
    tavg = (texp(1:end-1)+texp(2:end))/2;
%     Davg = (DexpSmooth(1:end-1)+DexpSmooth(2:end))/2;
%     ep = diff(DexpSmooth)./diff(texp).*-2./Davg;
    Davg = (Dexp(1:end-1)+Dexp(2:end))/2;
    ep = diff(Dexp)./diff(texp).*-2./Davg;
%     figure();
%     plot(tavg,ep)
    % linear interpolant for epsilon and a
    Fep = griddedInterpolant(tavg, ep);
    Fa = griddedInterpolant(texp, Dexp/2);
    tmax = 10;

    options = odeset('Events', @(t,y) stressBalanced(t,y,...
        chi, X, mu, rho, Nm, Ns, L, lambda, g, Fa, Fep), ...
        'RelTol', 1e-6,'AbsTol',1e-9);
    y0 = ones(Nm*Ns*2,1);
    
%     tic
    [t,y] = ode15s(@(t,y) odefun(t,y, Fep, Nm, Ns, L, lambda), ...
        [texp(1), tmax], y0, options);
%     toc
    
    t0 = t(end);
    a0 = Fa(t0);
    y0 = y(end,:);

end

function dydt = odefun(t, y, Fep, Nm, Ns, L, lambda)

    % manually
    count = 0;
    for mm=1:Nm
        for ss = 1:Ns
            count = count + 1;
            Az(mm,ss) = y(count);
            count = count + 1;
            Ar(mm,ss) = y(count);
        end
    end

    L2 = (L).^2;

    f = L2./(L2+3-sum(Az-2*Ar,2));

    ep = Fep(t);

    Azdot = 2*ep*Az - f./lambda.*(Az-1);
    Ardot = -ep*Ar - f./lambda.*(Ar-1);
    
    count = 0;
    for mm=1:Nm
        for ss = 1:Ns
            count = count + 1;
            dydt(count) = Azdot(mm,ss);
            count = count + 1;
            dydt(count) = Ardot(mm,ss);
        end
    end
    dydt = dydt';

end

function [value, isterminal, direction] = stressBalanced(t,y,...
        chi, X, mu, rho, Nm, Ns, L, lambda, g, Fa, Fep)

    count = 0;
    for mm=1:Nm
        for ss = 1:Ns
            count = count + 1;
            Az(mm,ss) = y(count);
            count = count + 1;
            Ar(mm,ss) = y(count);
        end
    end

    a = Fa(t);
    ep = Fep(t);

    L2 = (L).^2;

    f = L2./(L2+3-sum(Az-2*Ar,2));

    value = 1/8*ep^2*a^2*rho + sum(g.*f.*sum(Az-Ar,2)) ...
        + 3*mu*ep - (2*X-1)*chi/a;
    isterminal = 1;
    direction = 0;
end







