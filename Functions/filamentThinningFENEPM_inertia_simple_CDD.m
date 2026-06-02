function [a,Az,Ar,f,stress,ep,Wi,t,te,ye,ie] = ...
    filamentThinningFENEPM_inertia_simple_CDD(chi, X, mu, rho, ...
        Nm, Ns, Lm, Nk, zet_rat, lambda_ms, gm, a0, cutoff, tmax)

%     mvals = 1:Nm;
%     svals = 1:Ns;
%     indices = ((mvals-1)*Ns+svals');
%     indices = indices(:);

    Ntot = Nm*Ns;
    t0 = 0;

    y0(1) = a0;
    y0(2:2*Ntot+1) = 1;
    
    constants = [Nm, Ns, chi, mu, X, rho, a0];
    options = odeset('Events', @(t,y) aDecayed(t,y,a0,cutoff),...
        'RelTol', 1e-6,'AbsTol',1e-9);
    
    tic
    [t,y,te,ye,ie] = ode15s(@(t,y) odefun(t,y,constants,...
        Lm, lambda_ms, Nk, zet_rat, gm), ...
        [t0, tmax], y0, options);
    toc
    
    a = y(:,1);
    lent = length(t);
%     Az = reshape(y(:,2*indices),[lent,Nm,Ns]);
%     Ar = reshape(y(:,2*indices+1),[lent,Nm,Ns]);
    % manually
    count = 1;
    for mm=1:Nm
        for ss = 1:Ns
            count = count + 1;
            Az(:,mm,ss) = y(:,count);
            count = count + 1;
            Ar(:,mm,ss) = y(:,count);
        end
    end
    L2 = (Lm).^2;
%     f = sum((L2'*Ns)./((L2'*Ns)+3-Az-2*Ar),3);
%     f = (L2'*Ns)./((L2'*Ns)+3-sum(Az-2*Ar,3));
%      f = L2'./(L2'+3-sum(Az-2*Ar,3)/Ns);
     f = L2'./(L2'+3-sum(Az+2*Ar,3));
%     f = sum(L2'./(L2'+3-Az-2*Ar),3)/Ns;
    stress = gm'.*f.*sum(Az-Ar,3);
    ep = ones(size(t))*sqrt(8*chi/(a(1)^3*rho));
    if Ns>1
        for tt = 1:lent
            ep(tt) = fzero(@(ep) stressBalance(ep, a(tt), ...
                    squeeze(Az(tt,:,:)), squeeze(Ar(tt,:,:)),...
                    constants, Lm, gm), ep(tt));
        end
    else
        for tt = 1:lent
            ep(tt) = fzero(@(ep) stressBalance(ep, a(tt), ...
                    squeeze(Az(tt,:,:))', squeeze(Ar(tt,:,:))',...
                    constants, Lm, gm), ep(tt));
        end
    end
    Wi = lambda_ms(:,1).*ep';

end

function dydt = odefun(t, y, constants, L, lambda, Nk, zet_rat, g)
    Nm = constants(1);
    Ns = constants(2);
    chi = constants(3);
    mu = constants(4);
    X = constants(5);
    rho = constants(6);
    a0 = constants(7);

%     mvals = 1:Nm;
%     svals = 1:Ns;
%     indices = ((mvals-1)*Ns+svals');
%     indices = indices(:);

    a = y(1);
%     Azt = reshape(y(2*indices),[Ns,Nm])';
%     Art = reshape(y(2*indices+1),[Ns,Nm])';

    himodel = 2;
    % phieq = 0;
    d = 0;

    % manually
    count = 1;
    for mm=1:Nm
        for ss = 1:Ns
            count = count + 1;
            Az(mm,ss) = y(count);
            count = count + 1;
            Ar(mm,ss) = y(count);
            Q = sqrt((Az(mm,ss) + 2*Ar(mm,ss))/3);
            if ss == 1
                zet(mm,ss) = (Q-1)/L(mm)*(zet_rat-1) + 1;
            else 
                zet(mm,ss) = 1;
            end
        end
    end

    L2 = (L).^2;

%     f = sum((L2*Ns)./((L2*Ns)+3-Az-2*Ar),2);
%     f = (L2*Ns)./((L2*Ns)+3-sum(Az-2*Ar,2));
%     f = L2./(L2+3-sum(Az-2*Ar,2)/Ns);
    f = L2./(L2+3-sum(Az-2*Ar,2));
%     f = sum(L2./(L2+3-Az-2*Ar),2)/Ns;

    persistent epGuess
    if a==a0
        epGuess = sqrt(8*chi/(a^3*rho));
    end

    ep = fzero(@(ep) stressBalance(ep, a, Az, Ar, constants, L, g), epGuess);
%     e = (chi./a' - sum(g.*f.*sum(Az-Ar,2)))/(3*mu)
    
    epGuess = ep;

    adot = -0.5*ep*a;
    Azdot = 2*ep*Az - f./(lambda.*zet).*(Az-1);
    Ardot = -ep*Ar - f./(lambda.*zet).*(Ar-1);
    
    dydt(1) = adot;
    count = 1;
    for mm=1:Nm
        for ss = 1:Ns
            count = count + 1;
            dydt(count) = Azdot(mm,ss);
            count = count + 1;
            dydt(count) = Ardot(mm,ss);
        end
    end
%     dydt(2*indices) = Azdot(:);
%     dydt(2*indices+1) = Ardot(:);
    dydt = dydt';
%     dydt = [adot; Azdot(:); Ardot(:)];

end

function bal = stressBalance(ep, a, Az, Ar, constants, L, g)
    Nm = constants(1);
    Ns = constants(2);
    chi = constants(3);
    mu = constants(4);
    X = constants(5);
    rho = constants(6);
    a0 = constants(7);

    L2 = (L).^2;
    
%     f = (L2*Ns)./((L2*Ns)+3-sum(Az-2*Ar,2));
%     f = sum((L2*Ns)./((L2*Ns)+3-Az-2*Ar),2);
%     f = L2./(L2+3-sum(Az-2*Ar,2)/Ns);
    f = L2./(L2+3-sum(Az+2*Ar,2));
%     f = sum(L2./(L2+3-Az-2*Ar),2)/Ns;

    bal = 1/8*ep^2*a^2*rho + sum(g.*f.*sum(Az-Ar,2)) ...
        + 3*mu*ep - (2*X-1)*chi/a;
end


function [value, isterminal, direction] = aDecayed(t,y,a0,cutoff)
    a = y(1);
    value = a/a0-cutoff;
    isterminal = 1;
    direction = 0;
end











