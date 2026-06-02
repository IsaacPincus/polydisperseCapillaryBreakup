% clear variables;
% close all;
warning('off','MATLAB:nearlySingularMatrix');

datapath = "../Data";
addpath(datapath)
addpath("../Functions")

% spreadsheetName = 'ThinningCurves_PS7_PS16_andMixtures.xlsx';

% PS constants
alpha = 0.5;
K = 8e-8; % m^3/g, Mark-Houwink parameter, fitted
% alpha = 1;
% K = 5e-11; % m^3/g, Mark-Houwink parameter, fitted
density_solvent = 0.98*997e3; % g/m^3
vol_solvent = 1e6/density_solvent; % vol in m^3 of 1Mil g of solvent
eta_s = 0.059; % Pa.s, solvent viscosity
mu = eta_s;
% eta_s = 100; % Pa.s, solvent viscosity
NA = 6.02e23;  % 1/mol, Avagadro's constant
kB = 1.38e-23;  % J/k, Boltzmann constant
T = 298;        % K, temperature
l = 1.54e-10;   % m, PS bond length
cinf = 9.7;     % PS characteristic ratio
MW_monomer = 104; % g/mol, PS monomer weight
nu = 0.5;
chi = 3.1e-2; % N/m, DOP maybe 0.015 N/m?
rho = 1000; % kg/m^3, fluid density
alpha_C1 = 0.25;
C1 = 9/(4*alpha_C1^3);
rho = rho*C1*8;
tmax = 100; % s
cutoff = 1e-3;
hstar = 0.267;
hsk = 0.027;
% hsk = 0.0;
phieq = 0;
% Ctau = 1;
Uetatau = 2.36;
a0 = 6.3e-4;
X = 1;
Ns = 10;
PS7_conc = 200;
PS16_conc_vals = [20,100,200,500];
zeta0 = 0.5;
xi0 = 1;
R0 = 3e-3;
L0 = 1e-3;
V0_vals = linspace(0.1,1.5,10)*1e-3;
Lend = 1.3e-3;
Hend = Lend/L0;
Lam0 = L0/R0;
% linear
Hfunc = @(t) 1 + t;
Hderivfunc = @(t) 1;
% % exponential
% Hfunc = @(t) exp(t);
% Hderivfunc = @(t) exp(t);
colours_MWDist = parula(3);
colours_MWDist = colours_MWDist([1,2],:);
graphs_on = false;
% graphs_on = true;

% % PS constants
% alpha = 0.5;
% K = 11e-8; % m^3/g, Mark-Houwink parameter, fitted
% % alpha = 1;
% % K = 5e-11; % m^3/g, Mark-Houwink parameter, fitted
% density_solvent = 0.98*997e3; % g/m^3
% vol_solvent = 1e6/density_solvent; % vol in m^3 of 1Mil g of solvent
% eta_s = 0.056; % Pa.s, solvent viscosity
% % eta_s = 100; % Pa.s, solvent viscosity
% NA = 6.02e23;  % 1/mol, Avagadro's constant
% kB = 1.38e-23;  % J/k, Boltzmann constant
% T = 298;        % K, temperature
% l = 1.54e-10;   % m, PS bond length
% cinf = 9.7;     % PS characteristic ratio
% MW_monomer = 104; % g/mol, PS monomer weight
% nu = 0.5;
% chi = 3.1e-2; % N/m, DOP maybe 0.015 N/m?
% rho = 1000; % kg/m^3, fluid density
% alpha_C1 = 0.4;
% C1 = 9/(4*alpha_C1^3);
% rho = rho*C1*8;
% tmax = 100; % s
% cutoff = 1e-4;
% Ns = 10;
% PS7_conc = 200;
% PS16_conc_vals = [20,100,200,500];
% hstar = 0.2671;
% mu = eta_s;
% % a0 = 4.8e-4;
% X = 1;
% zeta0 = 0.5;
% xi0 = 1;
% R0 = 3e-3;
% L0 = 1e-3;
% V0_vals = linspace(0.1,1.5,10)*1e-3;
% Lend = 1.3e-3;
% Hend = Lend/L0;
% Lam0 = L0/R0;
% % linear
% Hfunc = @(t) 1 + t;
% Hderivfunc = @(t) 1;
% % % exponential
% % Hfunc = @(t) exp(t);
% % Hderivfunc = @(t) exp(t);
% colours_MWDist = parula(3);
% colours_MWDist = colours_MWDist([1,2],:);
% graphs_on = false;
% % graphs_on = true;

opts = detectImportOptions("PS7MWDist.csv");
dataPS7 = readmatrix("PS7MWDist.csv", opts);
PS7MW = dataPS7(:,1);
PS7NormMass = dataPS7(:,2);

opts = detectImportOptions("PS16MWDist.csv");
dataPS16 = readmatrix("PS16MWDist.csv", opts);
PS16MW = dataPS16(:,1);
PS16NormMass = dataPS16(:,2);

for jj = 1:length(PS16_conc_vals)
    PS16_conc = PS16_conc_vals(jj);
    for ii = 1:length(V0_vals)
        V0 = V0_vals(ii);
        tpull = fzero(@(t) Hfunc(t) - Hend, 1)*L0/V0;
    %     %% fit individual, two modes
        
        % project data onto same x-values
        allPS_data = [PS7MW;PS16MW];
        noPoints = 10000;
        x = linspace(min(allPS_data), max(allPS_data),noPoints);
        method = 'pchip';
        yPS7 = interp1(PS7MW,PS7NormMass,x, method, 0);
        yPS16 = interp1(PS16MW,PS16NormMass,x, method, 0);
        
        Nm = 1;
        mass_polymer = PS7_conc;
        
        weight_fraction_PS16 = 0/mass_polymer;
        y = yPS7*(1-weight_fraction_PS16)+yPS16*weight_fraction_PS16;
        density_mass_PS = mass_polymer/vol_solvent;
        maxMW = max(x);
        [nPS7, MPS7, mol_polymer_densityPS7] = ...
            splitMWDist(graphs_on, colours_MWDist,...
            x, y, noPoints, maxMW, Nm, nu, density_mass_PS);
        
        Nm = 1;
        mass_polymer = PS16_conc;
        
        weight_fraction_PS16 = PS16_conc/mass_polymer;
        y = yPS7*(1-weight_fraction_PS16)+yPS16*weight_fraction_PS16;
        density_mass_PS = mass_polymer/vol_solvent;
        maxMW = max(x);
        [nPS16, MPS16, mol_polymer_densityPS16] = ...
            splitMWDist(graphs_on, colours_MWDist,...
            x, y, noPoints, maxMW, Nm, nu, density_mass_PS);
        
        % figure(2);
        % hold on
        % copyobj(figure(2).Children.Children, figure(1).Children);
        % 
        % close 2
        % 
        % f1 = gcf;
        % objs = f1.Children.Children;
        % 
        % objs(1).Color = colours_MWDist(1,:);
        % objs(2).FaceColor = colours_MWDist(1,:);
        % % objs(3).Color = colours_MWDist(1,:);
        % % objs(3).LineWidth = 3;
        % 
        % objs(4).Color = colours_MWDist(2,:);
        % objs(5).FaceColor = colours_MWDist(2,:);
        % % objs(6).Color = colours_MWDist(2,:);
        % % objs(6).LineWidth = 3;
        
        Nm = 2;
        mol_polymer_density = mol_polymer_densityPS7 + mol_polymer_densityPS16;
        mol_frac_PS7 = mol_polymer_densityPS7/mol_polymer_density;
        n = [nPS16*(1-mol_frac_PS7), nPS7*mol_frac_PS7];
        M = [MPS16, MPS7];
        
        % objs(1).YData = [0,objs(4).YData(2)*n(1)/n(2)];
        
        % get densities etc
        density = mol_polymer_density; % mol/m^3, total mol density of polymer
        g = n'*density*kB*T*NA;
        
        Ltrue = 2*l*M/MW_monomer;
        R0_pol = sqrt(2*l^2.*M/MW_monomer*cinf);
        Rg = R0_pol/sqrt(6);
        Nk = Ltrue.^2./(2*cinf*l^2.*M/MW_monomer);
        Nk = Nk';
        L = sqrt(3*((Ltrue'./R0_pol').^2-1));
        
        jvals = 1:Ns;
        %     hstar = 0;
        sigma = -1.40*hstar^0.78;
        b = 1-1.66*hstar^0.78;
        intrinsic_visc = K*M.^alpha;
        lambda_m = M.*intrinsic_visc*eta_s/(NA*kB*T);
        %     lambda_m = M.^1.5.*intrinsic_visc*eta_s/(NA*kB*T);
        lambda_ms = lambda_m'./jvals.^(2+sigma);
        
    %     %% solve for stresses
        
        % firstly initial deformation
        tic
        options = odeset('RelTol', 1e-6,'AbsTol',1e-9);
        y0 = ones(Nm*Ns*2,1);
        
        sol = ode15s(@(t,y) odefun(t, y, Nm, Ns, L, lambda_ms,...
            Hfunc, Hderivfunc, zeta0, Lam0, V0, L0), ...
            [0, tpull], y0, options);
        toc
        
        points = 1000;
        [y,yp] = deval(sol,linspace(0,tpull,points));
        H = Hfunc(tpull*V0/L0);
        dH = Hderivfunc(tpull*V0/L0);
        z0 = ((zeta0-1)*zeta0)/(1-2*zeta0)^2;
        Rpull = (xi0*((H/z0 - 4)./(H/z0 - H*4)).^(3/(4*Lam0^2)))*R0
        Rpull*1000
        epend = (-3/(2*Lam0^2)*1./(H/(4*z0) - 1).*dH./H)*V0/L0
        ep0 = (-3/(2*Lam0^2)*1./(Hfunc(0)/(4*z0) - 1).*Hderivfunc(0)./Hfunc(0))*V0/L0
        
        % manually
        Az = zeros(Nm,Ns,points);
        Ar = zeros(Nm,Ns,points);
        count = 0;
        for mm=1:Nm
            for ss = 1:Ns
                count = count + 1;
                Az(mm,ss,:) = y(count,:);
                count = count + 1;
                Ar(mm,ss,:) = y(count,:);
            end
        end
        
        L2 = (L).^2;
        f = L2./(L2+3-sum(Az+2*Ar,2));
        stressPolymerPull = sum(g.*f.*sum(Az-Ar,2));
        stressPolymerPull(end)
        stressCapPull = - (2*X-1)*chi/Rpull
        
    %     %%
        
        Ntot = Nm*Ns
        a0 = Rpull;
        epInitial = epend;
        y0(1) = a0;
        y0(2:2*Ntot+1) = y(:,end);
        
        [a,Az,Ar,f,stress,ep,Wi,t,te,ye,ie] = ...
            filamentThinningFENEPM_inertia_initial(chi, X, mu, rho, ...
            Nm, Ns, L, lambda_ms, g, a0, cutoff, tmax, y0, epInitial);
        
        local_ep_maxes = find(islocalmax(ep)==1);
        t_max_ep_all = t(islocalmax(ep)==1);
    %     t_c = t_max_ep_all(1);
        t_c = 0;
        
    %     ep_min_array = islocalmin(ep)==1;
    %     ep_min_all = ep(ep_min_array&(t>t_c));
    %     ep_min = min(ep_min_all);
    %     t_ep_min = min(t(ep_min_array&(t>t_c)));
        [ep_min, index_ep_min] = min(ep);
        t_ep_min = t(index_ep_min);
        
        tau_ep_min(ii,jj) = 2/(3*ep_min);
        
        % indices_EC = find((ep<1.2*ep_min)&(t>t_c));
        indices_EC = find((ep<1.2*ep_min)&(t>t_c));
        [xxx,yyy] = prepareCurveData(t(indices_EC)-t_c,a(indices_EC));
        fitobject = fit(xxx, yyy, 'exp1');
        tau_EC_fit = -1./(3*fitobject.b);
        tau_EC(ii,jj) = tau_EC_fit;
        
        %     t_c(ii) = 0;
        
        aModel = a;
        epModel = ep;
        tModel = t;
        stressModel = stress;
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting tau as a function of V0
figure();
hold on
plot(V0_vals, tau_EC)
% plot(V0_vals, tau_ep_min)
axes1 = gca;
axes1.XScale = 'log';
xlabel('$V_0$ [mm/s]')
ylabel('$\lambda$ [s]')
ylim([0,0.16])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plotting fits
% fig = figure();
% fig.Position = [488,152.2,810,610];
% hold on
% axes1 = gca;
% axes1.YScale = 'log';
% tm = t-t_c;
% plot(fitobject)
% plot(tm, a, 'k');
% plot(tm(indices_EC), a(indices_EC), 'or');
% legend({'Exponential Fit', 'Simulation', 'Points used for fit'})
% xlabel('$t-t_c$ [s]');
% ylabel('$a$ [m]');
% xlim([-0.015,1.2]);
% ylim([1e-6,3e-4]);
% 
% fig = figure();
% fig.Position = [488,152.2,810,610];
% hold on
% axes1 = gca;
% axes1.YScale = 'log';
% tm = t-t_c;
% plot(tm, ep, 'k');
% plot(t_ep_min-t_c, ep_min, 'or');
% plot(tm, ones(size(ep))*ep_min, 'b--', 'linewidth', 0.5);
% plot(tm, ones(size(ep))*ep_min*1.2, 'b--', 'linewidth', 0.5);
% legend({'Simulation', 'Min $\dot{\epsilon}$'})
% xlabel('$t-t_c$ [s]');
% ylabel('$\dot{\epsilon}$ [s$^{-1}$]');
% xlim([-0.015,1.2]);
% ylim([4,1000]);

%% plotting individual for last run

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % epsilon
% figure();
% hold on
% plot(tModel-t_c, epModel)
% 
% xlabel('$t-t_c$ [s]')
% ylabel('$\dot{\epsilon}$')
% % legend(string(1:5))
% % colororder(colours_MWDist);
% 
% % axes1 = gca;
% % axes1.YScale = 'log';
% % axes1.XScale = 'log';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % contributions to stress
% fig = figure();
% fig.Position = [488,152.2,810,610];
% hold on
% 
% % capillary stress
% plot(tModel-t_c, chi./aModel, 'k', 'linewidth', 2);
% 
% xlabel('$t-t_c$ [s]')
% % ylabel('$g_i f_i (A_{z,i} - A_{r,i})$')
% ylabel('stress [Pa]')
% 
% axes1 = gca;
% axes1.YScale = 'log';
% % axes1.XScale = 'log';
% 
% % polymer stresses
% for pp = 1:Nm
%     plot(tModel-t_c, stress(:,pp), '-', 'Color',colours_MWDist(pp,:));
%     heightIndex = find(tModel-t_c>2*lambda_ms(pp),1);
%     plot(2*[lambda_ms(pp),lambda_ms(pp)], ...
%           [1e-5,stress(heightIndex,pp)], ':',...
%           'HandleVisibility','off',...
%           'Color',colours_MWDist(pp,:))
% end
% 
% % intertial and viscous stresses
% plot(tModel-t_c, 3*eta_s*ep + 1/8*rho*ep.^2.*aModel.^2, ...
%     'k--', 'linewidth', 1.5);
% 
% l1 = legend({'Capillary Stress',...
%     sprintf('%0.3g MDa Polymer Stress', M(1)/1e6),...
%     sprintf('%0.3g MDa Polymer Stress', M(2)/1e6),...
%     'Viscous and Inertial Stresses'});
% 
% xlim([tModel(1)-t_c,tModel(end)-t_c])
% ylim([0.5,1e5])
% 
% copyobj(fig.Children(2), fig);
% axes2 = fig.Children(3);
% axes2.Position = [0.42,0.39,0.218,0.23];
% axes2.XTick = [];
% axes2.YTick = [];
% xlabel('')
% ylabel('')
% xlim([tModel(1)-t_c,0.01])
% ylim([30,500])
% fig.Children = flip(fig.Children);
% 
% l1.Position = [0.162,0.719,0.388,0.174];
% 
% % Create arrow
% annotation(fig,'arrow',[0.41600790513834,0.141304347826087],...
%     [0.55570117955439,0.486238532110092]);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Az/L^2
% figure();
% hold on
% 
% yyaxis right
% plot(tModel-t_c, ep, 'k-', 'linewidth', 2)
% ylabel('$\dot{\epsilon}$')
% 
% yyaxis left
% scaledAz = sum(Az,3)./L'.^2;
% p1 = plot(tModel-t_c, scaledAz, '-');
% 
% xlabel('$t-t_c$ [s]')
% ylabel('$A_{z,i}/L_i^2$')
% % legend(string(1:5))
% colororder(colours_MWDist);
% 
% axes1 = gca;
% axes1.YScale = 'log';
% % axes1.XScale = 'log';
% axes1.YAxis(1).Color = 'k';
% axes1.YAxis(2).Color = 'k';
% 
% % add labels
% for pp = 1:length(p1)
%     offset = 0.65+(pp-length(p1)/2)^2/(0.6*length(p1)^2);
% %     offset = 1;
%     label(p1(pp), sprintf('%0.3g MDa', M(pp)/1e6),...
%         'location','center', 'offset', offset)
%     heightIndex = find(tModel-t_c>2*lambda_ms(pp),1);
%     plot(2*[lambda_ms(pp),lambda_ms(pp)], ...
%       [1e-5,scaledAz(heightIndex,pp)], ':',...
%         'HandleVisibility','off')
% end
% 
% % Az
% figure();
% hold on
% plot(tModel-t_c, sum(Az, 3)/Ns)
% 
% xlabel('$t-t_c$ [s]')
% ylabel('$A_{z,i}$')
% % legend(string(1:5))
% colororder(colours_MWDist);
% 
% axes1 = gca;
% axes1.YScale = 'log';
% % axes1.XScale = 'log';
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % extension rate
% fig = figure();
% fig.Position = [488,152.2,810,610];
% hold on
% 
% for pp = 1:Nm
%     plot(tModel-t_c, epModel*lambda_m(pp),...
%         '-', 'Color',colours_MWDist(pp,:), ...
%         'DisplayName',sprintf('%0.3g MDa', M(pp)/1e6));
% end
% 
% xlabel('$t-t_c$ [s]')
% ylabel('$W\!i$')
% legend('Location','northeast')
% 
% axes1 = gca;
% axes1.YScale = 'log';
% % axes1.XScale = 'log';
% 
% xlim([tModel(1)-t_c,tModel(end)-t_c])
% ylim([1e-1,1e2])
% 
% plot(tModel-t_c, ones(size(tModel))*0.5, 'k--', ...
%     'DisplayName','$Wi = 0.5$ (coil-stretch)')
% 
% copyobj(fig.Children(2), fig);
% axes2 = fig.Children(1);
% axes2.Position = [0.254,0.586,0.218,0.230];
% % axes2.XTick = [];
% % axes2.YTick = [];
% xlabel(axes2, '')
% ylabel(axes2, '')
% xlim(axes2, [tModel(1)-t_c,0.01])
% ylim(axes2, [1e0,1e2])
% 
% % fig.Children = flip(fig.Children);
% % Create arrow
% annotation(fig,'arrow',[0.253086419753086,0.144444444444444],...
%     [0.783262295081967,0.768524590163934]);


% %% extension and f-factor
% % Az
% figure();
% hold on
% plot(tModel-t_c, sum(Az, 3)/Ns)
% 
% xlabel('$t-t_c$ [s]')
% ylabel('$A_{z,i}$')
% % legend(string(1:5))
% colororder(colours_MWDist);
% 
% axes1 = gca;
% axes1.YScale = 'log';
% % axes1.XScale = 'log';
% 
% figure();
% hold on
% plot(tModel-t_c, f)
% % axes1 = gca;
% % axes1.YScale = 'log';

function dydt = odefun(t, y, Nm, Ns, L, lambda,...
    Hfunc, Hderivfunc, zeta0, Lam0, V0, L0)

    % general analytical solution for filament thinning prior to
    % instability
    H = Hfunc(t*V0/L0);
    dH = Hderivfunc(t*V0/L0);
    z0 = ((zeta0-1)*zeta0)/(1-2*zeta0)^2;
    ep = (-3/(2*Lam0^2)*1./(H/(4*z0) - 1).*dH./H)*V0/L0;

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
    f = L2./(L2+3-sum(Az+2*Ar,2));

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

