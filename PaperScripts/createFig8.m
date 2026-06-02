clear variables;
close all;
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
Nm = 80;
Ns = 1;
hstar = 0.267;
hsk = 0.027;
% hsk = 0.0;
phieq = 0;
% Ctau = 1;
Uetatau = 2.36;
a0 = 6.3e-4;
X = 1;
PS7_conc_vals = [200,500,1000];
PS16_conc_vals = linspace(0,300,50);
colours_MWDist = parula(Nm);
colours_PSvals = winter(length(PS16_conc_vals));
graphs_on = false;

opts = detectImportOptions("PS7MWDist.csv");
dataPS7 = readmatrix("PS7MWDist.csv", opts);
PS7MW = dataPS7(:,1);
PS7NormMass = dataPS7(:,2);

opts = detectImportOptions("PS16MWDist.csv");
dataPS16 = readmatrix("PS16MWDist.csv", opts);
PS16MW = dataPS16(:,1);
PS16NormMass = dataPS16(:,2);

%% fit together, multiple modes, Fig8a

for jj = 1:length(PS7_conc_vals)
    PS7_conc = PS7_conc_vals(jj);
    for ii = 1:length(PS16_conc_vals)
        PS16_conc = PS16_conc_vals(ii);
        mass_polymer = PS16_conc + PS7_conc;
        
        % project data onto same x-values
        allPS_data = [PS7MW;PS16MW];
        noPoints = 10000;
        x = linspace(min(allPS_data), max(allPS_data),noPoints);
        method = 'pchip';
        yPS7 = interp1(PS7MW,PS7NormMass,x, method, 0);
        yPS16 = interp1(PS16MW,PS16NormMass,x, method, 0);
        
        weight_fraction_PS16 = PS16_conc/mass_polymer;
        y = yPS7*(1-weight_fraction_PS16)+yPS16*weight_fraction_PS16;
        density_mass_PS = mass_polymer/vol_solvent;
        maxMW = max(x);
        [n, M, mol_polymer_density] = ...
            splitMWDist(graphs_on, colours_MWDist,...
            x, y, noPoints, maxMW, Nm, nu, density_mass_PS);
        
    %     %% solve for stresses
        % get densities etc
        density = mol_polymer_density; % mol/m^3, total mol density of polymer
        g = n'*density*kB*T*NA;
        
        Ltrue = 2*l*M/MW_monomer;
        R0 = sqrt(2*l^2.*M/MW_monomer*cinf);
        Rg = R0/sqrt(6);
        Nk = Ltrue.^2./(2*cinf*l^2.*M/MW_monomer);
        Nk = Nk';
        L = Ltrue'./R0';
        
        jvals = 1:Ns;
        %     hstar = 0;
        sigma = -1.40*hstar^0.78;
        b = 1-1.66*hstar^0.78;
        intrinsic_visc = K*M.^alpha;
        lambda_m = M.*intrinsic_visc*eta_s/(NA*kB*T);
        %     lambda_m = M.^1.5.*intrinsic_visc*eta_s/(NA*kB*T);
        lambda_ms = lambda_m'./jvals.^(2+sigma);
        
        [a,~,~,~,stress,ep,Wi,t,~,~,~] = ...
            filamentThinningFENEPM_inertia_CDD(chi, X, mu, rho, ...
            Nm, Ns, L, Nk, hsk, lambda_ms, g, a0, cutoff, tmax, phieq);
        
        local_ep_maxes = find(islocalmax(ep)==1);
        t_max_ep_all = t(islocalmax(ep)==1);
        t_c{ii,jj} = t_max_ep_all(1);
        
        ep_min_array = islocalmin(ep)==1;
        indices_ep_min = find(ep_min_array&(t>t_c{ii,jj}));
        ep_min_all = ep(indices_ep_min);
        [ep_min{ii,jj}, index] = min(ep_min_all);
        ep_min_index{ii,jj} = indices_ep_min(index);
        
        tau_ep_min{ii,jj} = 2/(3*ep_min{ii,jj});
        
        indices_EC = find((ep<2*ep_min{ii,jj})&(t>t_c{ii,jj}));
        [xxx,yyy] = prepareCurveData(t(indices_EC),a(indices_EC));
        fitobject = fit(xxx, yyy, 'exp1');
        tau_EC_fit = -1./(3*fitobject.b);
        tau_EC{ii,jj} = tau_EC_fit;
        
        aModel{ii,jj} = a;
        epModel{ii,jj} = ep;
        tModel{ii,jj} = t;
        stressModel{ii,jj} = stress;
        
        nModel{ii,jj} = n;
        MModel{ii,jj} = M;
        lambdaModel{ii,jj} = lambda_m;
        gModel{ii,jj} = g;
        LModel{ii,jj} = L;
    
    end
end

%% plotting individual for last run

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % contribution to stress at minimum epsilon
% 
% fig = figure();
% fig.Position = [488,152.2,810,610];
% hold on
% 
% for ii = 1:length(PS16_conc_vals)
%     stress_min = stressModel{ii,jj}(ep_min_index{ii,jj},:);
%     Mm = MModel{ii,jj}*1e-6;
%     plot(Mm, -stress_min./trapz(Mm,stress_min), '-', ...
%         'DisplayName','$\sigma_p$-weighted',...
%         'Color',colours_PSvals(ii,:))
%     plot(Mm, -nModel{ii,jj}./trapz(Mm,nModel{ii,jj}), ':', ...
%         'DisplayName','$n$-weighted',...
%         'Color',colours_PSvals(ii,:))
% end
% % legend
% xlabel('MW [10$^{6}$ g/mol]')
% ylabel('pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Same but with times rather than MW

% fig = figure();
% fig.Position = [488,152.2,810,610];
% hold on

for jj = 1:length(PS7_conc_vals)
    for ii = 1:length(PS16_conc_vals)
        stress_min = stressModel{ii,jj}(ep_min_index{ii,jj},:);
%         xs = lambdaModel{ii,jj};
        xs = MModel{ii,jj}*1e-6;
        ys = -stress_min./trapz(xs,stress_min);
%         plot(xs, ys, '-', 'DisplayName','$\sigma_p$-weighted',...
%             'Color',colours_PSvals(ii,:))
%         plot(lambdaModel{ii,jj},...
%             -nModel{ii,jj}./trapz(lambdaModel{ii,jj},nModel{ii,jj}),...
%             ':', 'DisplayName','$n$-weighted',...
%             'Color',colours_PSvals(ii,:))
        meanStressP(ii,jj) = -trapz(xs,xs.*ys);
        [val, index] = max(ys);
%         modex(ii,jj) = xs(index);
        
        stress_ana = 2*gModel{ii,jj}.*lambdaModel{ii,jj}'...
            .*LModel{ii,jj}.^2.*epModel{ii,jj}(ep_min_index{ii,jj});
        ys = -stress_ana./trapz(xs,stress_ana);
        meanStressAna(ii,jj) = -trapz(xs,xs.*ys');

        meanNumber(ii,jj) = sum(nModel{ii,jj})/sum(nModel{ii,jj}./xs);

        tau_EC_val(ii,jj) = tau_EC{ii,jj};

    %     plot([meanx(ii),meanx(ii)],[0,10], 'r-', 'LineWidth', 2,...
    %         'DisplayName','$\dot{\epsilon}_\mathrm{min}$',...
    %         'Color',colours_PSvals(ii,:))
    
    end
end

% legend
% xlabel('$\lambda_{MW}$ [s]')
% ylabel('pdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot mean value against PS16 concentration
fig = figure();
fig.Position = [488,152.2,810,610];
hold on

plot(PS16_conc_vals,meanStressP, '-', 'DisplayName','Polymer')
plot(PS16_conc_vals,meanStressAna, ':', 'DisplayName','Analytical')
plot(PS16_conc_vals,meanNumber, '--', 'DisplayName','Number')
% plot(PS16_conc_vals,modex, 'k:', 'DisplayName','mode')
colororder(flip(copper(3)))

xlabel('$c_\mathrm{PS16}$ [ppm]')
ylabel('$\lambda$ [s]')
ylabel('$M$ [MDa]')

ylim([5.5,17])

legend([" "," "," "," "," "," ",...
    "$\qquad$200","$\qquad$500","$\qquad$1000"],...
    'Position',[0.4327 0.29017 0.40472 0.1499],...
    'NumColumns',3)

annotation('textbox',...
    [0.42888 0.43309 0.46 0.066229],...
    'String',['$\sigma_p$ $\qquad$ $\sigma_\mathrm{ana}$  $\qquad$  $n$' ...
    '  $\quad$   $c_\mathrm{PS7}$ [ppm]'],...
    'FitBoxToText','off',...
    'BackgroundColor',[1,1,1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Same but with relaxation times and tau_EC
fig = figure();
fig.Position = [488,152.2,810,610];
hold on

plot(PS16_conc_vals,K*(meanStressP*1e6).^(1 + alpha)*eta_s/(NA*kB*T), ...
    '-', 'DisplayName','Polymer')
plot(PS16_conc_vals,K*(meanStressAna*1e6).^(1 + alpha)*eta_s/(NA*kB*T), ...
    ':', 'DisplayName','Analytical')
plot(PS16_conc_vals,K*(meanNumber*1e6).^(1 + alpha)*eta_s/(NA*kB*T), ...
    '--', 'DisplayName','Number')
plot(PS16_conc_vals, tau_EC_val, 'o', 'DisplayName','$\tau_\mathrm{EC}$')
% plot(PS16_conc_vals,modex, 'k:', 'DisplayName','mode')
colororder(flip(copper(3)))

xlabel('$c_\mathrm{PS16}$ [ppm]')
ylabel('$\lambda$ [s]')
% ylabel('$M$ [MDa]')

% ylim([5.5,15])

legend([" "," "," "," "," "," ",...
    "$\qquad$200","$\qquad$500","$\qquad$1000"],...
    'Position',[0.4327 0.29017 0.40472 0.1499],...
    'NumColumns',3)

annotation('textbox',...
    [0.42888 0.43309 0.46 0.066229],...
    'String',['$\sigma_p$ $\qquad$ $\sigma_\mathrm{ana}$  $\qquad$  $n$' ...
    '  $\quad$   $c_\mathrm{PS7}$ [ppm]'],...
    'FitBoxToText','off',...
    'BackgroundColor',[1,1,1]);

%% fit together, multiple modes, Fig8b

% change concentrations to what we want
PS7_conc = 1000;
PS16_conc = 100;
graphs_on = false;

mass_polymer = PS16_conc + PS7_conc;

% project data onto same x-values
allPS_data = [PS7MW;PS16MW];
noPoints = 10000;
x = linspace(min(allPS_data), max(allPS_data),noPoints);
method = 'pchip';
yPS7 = interp1(PS7MW,PS7NormMass,x, method, 0);
yPS16 = interp1(PS16MW,PS16NormMass,x, method, 0);

weight_fraction_PS16 = PS16_conc/mass_polymer;
y = yPS7*(1-weight_fraction_PS16)+yPS16*weight_fraction_PS16;
density_mass_PS = mass_polymer/vol_solvent;
maxMW = max(x);
[n, M, mol_polymer_density] = ...
    splitMWDist(graphs_on, colours_MWDist,...
    x, y, noPoints, maxMW, Nm, nu, density_mass_PS);

ML = (trapz(x, y.*x.^(1+nu))/trapz(x, y))^(1/(1+nu));


% %% solve for stresses
% get densities etc
density = mol_polymer_density; % mol/m^3, total mol density of polymer
g = n'*density*kB*T*NA;

Ltrue = 2*l*M/MW_monomer;
R0 = sqrt(2*l^2.*M/MW_monomer*cinf);
Rg = R0/sqrt(6);
Nk = Ltrue.^2./(2*cinf*l^2.*M/MW_monomer);
Nk = Nk';
L = Ltrue'./R0';

jvals = 1:Ns;
%     hstar = 0;
sigma = -1.40*hstar^0.78;
b = 1-1.66*hstar^0.78;
intrinsic_visc = K*M.^alpha;
lambda_m = M.*intrinsic_visc*eta_s/(NA*kB*T);
%     lambda_m = M.^1.5.*intrinsic_visc*eta_s/(NA*kB*T);
lambda_ms = lambda_m'./jvals.^(2+sigma);

[a,~,~,~,stress,ep,Wi,t,~,~,~] = ...
    filamentThinningFENEPM_inertia_CDD(chi, X, mu, rho, ...
    Nm, Ns, L, Nk, hsk, lambda_ms, g, a0, cutoff, tmax, phieq);

local_ep_maxes = find(islocalmax(ep)==1);
t_max_ep_all = t(islocalmax(ep)==1);
t_c = t_max_ep_all(1);

ep_min_array = islocalmin(ep)==1;
indices_ep_min = find(ep_min_array&(t>t_c));
ep_min_all = ep(indices_ep_min);
[ep_min, index] = min(ep_min_all);
ep_min_index = indices_ep_min(index);

tau_ep_min = 2/(3*ep_min);

indices_EC = find((ep<2*ep_min)&(t>t_c));
[xxx,yyy] = prepareCurveData(t(indices_EC),a(indices_EC));
fitobject = fit(xxx, yyy, 'exp1');
tau_EC_fit = -1./(3*fitobject.b);
tau_EC = tau_EC_fit;

%     t_c(ii) = 0;

aModel = a;
epModel = ep;
tModel = t;
stressModel = stress;

% %% plotting individual for last run

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% contribution to stress at minimum epsilon
stress_min_ana = 2*g.*lambda_m'.*L.^2.*ep_min;
stress_min = stress(ep_min_index,:);
xte = tModel(ep_min_index)-t_c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% overall figure with molecular weight
X1 = M*1e-6;
Ym1a = -stress_min_ana./trapz(M,stress_min_ana);
Ym1b = -stress_min./trapz(M,stress_min);
Ym1c = -n./trapz(M,n);
YMatrix1 = [Ym1a'; Ym1b; Ym1c];

X2 = tModel-t_c;
stress_cap = chi./aModel;
stress_other = 3*eta_s*ep + 1/8*rho*ep.^2.*aModel.^2;
YMatrix2 = [stress_cap, stress_other, stress];
X3 = [xte,xte];
Y1 = [1e-5,1e5];

createfigureMW(X1, YMatrix1, X2, YMatrix2, X3, Y1)



function createfigureMW(X1, YMatrix1, X2, YMatrix2, X3, Y1)
%CREATEFIGURE(X1, YMatrix1, X2, YMatrix2, X3, Y1)
%  X1:  vector of plot x data
%  YMATRIX1:  matrix of plot y data
%  X2:  vector of semilogy x data
%  YMATRIX2:  matrix of semilogy y data
%  X3:  vector of semilogy x data
%  Y1:  vector of semilogy y data

%  Auto-generated by MATLAB on 31-Oct-2024 10:40:38

% Create figure
figure1 = figure;
figure1.Position =[488,152.2,810,610];

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
ylim(axes1,[0, 3.5e-7]);

% Create multiple line objects using matrix input to plot
plot1 = plot(X1,YMatrix1,'Parent',axes1);
set(plot1(1),'DisplayName','$\sigma_{p,\mathrm{ana}}$-weighted',...
    'Color',[0 0 0], 'LineStyle',':');
set(plot1(2),'DisplayName','$\sigma_p$-weighted','Color',[0 0 0],...
    'LineStyle','-');
set(plot1(3),'DisplayName','$n$-weighted','Color',[0 0 0], ...
    'LineStyle','--');

% Create ylabel
ylabel('pdf');

% Create xlabel
xlabel('$M$ [MDa]');

hold(axes1,'off');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.145169485784945,0.756049340785143,0.253417204319446,0.149541631316831]);

% Create axes
axes2 = axes('Parent',figure1,...
    'Position',[0.531715463816916,0.597808865966411,0.281561264822132,0.298693750105661]);
hold(axes2,'on');
axes2.FontSize = 13;

% Create multiple line objects using matrix input to semilogy
semilogy1 = semilogy(X2,YMatrix2);
set(semilogy1(1),'DisplayName','Capillary','LineWidth',2,'Color',[0 0 0]);
set(semilogy1(2),'DisplayName','Viscous \& Inertial','LineStyle','--',...
    'Color',[0 0 0]);
set(semilogy1(3),'Color',[0.2422 0.1504 0.6603]);
set(semilogy1(4),...
    'Color',[0.248832911392405 0.16163417721519 0.698606329113924]);
set(semilogy1(5),...
    'Color',[0.254820253164557 0.17506835443038 0.733475949367089]);
set(semilogy1(6),...
    'Color',[0.260562025316456 0.18806582278481 0.768450632911392]);
set(semilogy1(7),...
    'Color',[0.265958227848101 0.200754430379747 0.803360759493671]);
set(semilogy1(8),...
    'Color',[0.270581012658228 0.214440506329114 0.835894936708861]);
set(semilogy1(9),...
    'Color',[0.274330379746835 0.229835443037975 0.864136708860759]);
set(semilogy1(10),...
    'Color',[0.277075949367089 0.246712658227848 0.888045569620253]);
set(semilogy1(11),...
    'Color',[0.279293670886076 0.264325316455696 0.908354430379747]);
set(semilogy1(12),...
    'Color',[0.280615189873418 0.282183544303797 0.925753164556962]);
set(semilogy1(13),...
    'Color',[0.281327848101266 0.30003164556962 0.940925316455696]);
set(semilogy1(14),...
    'Color',[0.281198734177215 0.31773417721519 0.954374683544304]);
set(semilogy1(15),...
    'Color',[0.280406329113924 0.335264556962025 0.966096202531646]);
set(semilogy1(16),...
    'Color',[0.278426582278481 0.352694936708861 0.976186075949367]);
set(semilogy1(17),...
    'Color',[0.275153164556962 0.370344303797468 0.984417721518987]);
set(semilogy1(18),...
    'Color',[0.270206329113924 0.388339240506329 0.990384810126582]);
set(semilogy1(19),...
    'Color',[0.263021518987342 0.406744303797468 0.994210126582279]);
set(semilogy1(20),...
    'Color',[0.252193670886076 0.42536582278481 0.997286075949367]);
set(semilogy1(21),...
    'Color',[0.236312658227848 0.444307594936709 0.999589873417722]);
set(semilogy1(22),...
    'Color',[0.216892405063291 0.464040506329114 0.996141772151899]);
set(semilogy1(23),...
    'Color',[0.197481012658228 0.483530379746835 0.989541772151899]);
set(semilogy1(24),...
    'Color',[0.185173417721519 0.501816455696202 0.982437974683544]);
set(semilogy1(25),...
    'Color',[0.180586075949367 0.519167088607595 0.973763291139241]);
set(semilogy1(26),...
    'Color',[0.177679746835443 0.536150632911392 0.963162025316456]);
set(semilogy1(27),...
    'Color',[0.175978481012658 0.552741772151899 0.949586075949367]);
set(semilogy1(28),...
    'Color',[0.169511392405063 0.568981012658228 0.936793670886076]);
set(semilogy1(29),...
    'Color',[0.157581012658228 0.585027848101266 0.92566582278481]);
set(semilogy1(30),...
    'Color',[0.148941772151899 0.600398734177215 0.914183544303797]);
set(semilogy1(31),...
    'Color',[0.143954430379747 0.615208860759494 0.903936708860759]);
set(semilogy1(32),...
    'Color',[0.136537974683544 0.62993417721519 0.896163291139241]);
set(semilogy1(33),...
    'Color',[0.125875949367089 0.644559493670886 0.889078481012658]);
set(semilogy1(34),...
    'Color',[0.114991139240506 0.658672151898734 0.880184810126582]);
set(semilogy1(35),...
    'Color',[0.103664556962025 0.672164556962025 0.868277215189873]);
set(semilogy1(36),...
    'Color',[0.0883379746835443 0.684720253164557 0.853656962025316]);
set(semilogy1(37),...
    'Color',[0.0649518987341772 0.696339240506329 0.836867088607595]);
set(semilogy1(38),...
    'Color',[0.0330658227848101 0.707018987341772 0.818449367088608]);
set(semilogy1(39),...
    'Color',[0.0081886075949367 0.716887341772152 0.799044303797468]);
set(semilogy1(40),...
    'Color',[0.0012873417721519 0.725962025316456 0.778788607594937]);
set(semilogy1(41),...
    'Color',[0.0138759493670886 0.734345569620253 0.758021518987342]);
set(semilogy1(42),...
    'Color',[0.0488620253164556 0.742037974683544 0.736663291139241]);
set(semilogy1(43),...
    'Color',[0.091312658227848 0.749139240506329 0.714825316455696]);
set(semilogy1(44),...
    'Color',[0.128686075949367 0.755883544303797 0.692641772151899]);
set(semilogy1(45),...
    'Color',[0.157586075949367 0.762596202531646 0.670112658227848]);
set(semilogy1(46),...
    'Color',[0.178725316455696 0.769474683544304 0.646878481012658]);
set(semilogy1(47),'Color',[0.195524050632911 0.776555696202532 0.6223]);
set(semilogy1(48),...
    'Color',[0.213318987341772 0.78333164556962 0.596173417721519]);
set(semilogy1(49),...
    'Color',[0.235792405063291 0.78956582278481 0.568218987341772]);
set(semilogy1(50),...
    'Color',[0.266413924050633 0.794692405063291 0.53873670886076]);
set(semilogy1(51),...
    'Color',[0.304421518987342 0.798436708860759 0.508007594936709]);
set(semilogy1(52),...
    'Color',[0.344440506329114 0.801015189873418 0.475672151898734]);
set(semilogy1(53),...
    'Color',[0.384444303797468 0.802717721518987 0.441083544303797]);
set(semilogy1(54),...
    'Color',[0.426959493670886 0.802875949367089 0.405439240506329]);
set(semilogy1(55),...
    'Color',[0.47309746835443 0.801206329113924 0.370710126582279]);
set(semilogy1(56),...
    'Color',[0.519778481012658 0.798108860759494 0.336587341772152]);
set(semilogy1(57),...
    'Color',[0.565122784810126 0.793874683544304 0.301688607594937]);
set(semilogy1(58),...
    'Color',[0.609683544303797 0.788643037974683 0.267070886075949]);
set(semilogy1(59),...
    'Color',[0.653477215189874 0.782305063291139 0.235160759493671]);
set(semilogy1(60),...
    'Color',[0.696135443037975 0.775029113924051 0.207601265822785]);
set(semilogy1(61),...
    'Color',[0.73706835443038 0.767362025316456 0.183693670886076]);
set(semilogy1(62),...
    'Color',[0.776272151898734 0.759292405063291 0.163584810126582]);
set(semilogy1(63),...
    'Color',[0.813513924050633 0.751289873417721 0.153862025316456]);
set(semilogy1(64),...
    'Color',[0.848426582278481 0.743722784810127 0.156217721518987]);
set(semilogy1(65),...
    'Color',[0.881653164556962 0.736959493670886 0.165569620253165]);
set(semilogy1(66),...
    'Color',[0.912525316455696 0.731639240506329 0.184158227848101]);
set(semilogy1(67),...
    'Color',[0.940532911392405 0.728708860759494 0.211340506329114]);
set(semilogy1(68),...
    'Color',[0.967324050632911 0.729067088607595 0.236141772151899]);
set(semilogy1(69),...
    'Color',[0.990174683544304 0.736659493670886 0.242767088607595]);
set(semilogy1(70),...
    'Color',[0.997126582278481 0.752788607594937 0.229783544303797]);
set(semilogy1(71),...
    'Color',[0.996751898734177 0.771113924050633 0.215977215189873]);
set(semilogy1(72),...
    'Color',[0.995122784810127 0.789856962025316 0.202341772151899]);
set(semilogy1(73),...
    'Color',[0.990686075949367 0.809191139240506 0.190767088607595]);
set(semilogy1(74),...
    'Color',[0.983003797468355 0.829098734177215 0.18116835443038]);
set(semilogy1(75),...
    'Color',[0.97366582278481 0.849211392405063 0.171344303797468]);
set(semilogy1(76),...
    'Color',[0.96583417721519 0.869324050632912 0.161827848101266]);
set(semilogy1(77),...
    'Color',[0.960953164556962 0.889336708860759 0.153548101265823]);
set(semilogy1(78),...
    'Color',[0.959517721518987 0.908922784810127 0.144716455696203]);
set(semilogy1(79),...
    'Color',[0.961116455696203 0.928035443037975 0.13313417721519]);
set(semilogy1(80),...
    'Color',[0.965016455696202 0.94680253164557 0.118987341772152]);
set(semilogy1(81),...
    'Color',[0.970667088607595 0.965378481012658 0.101467088607595]);
set(semilogy1(82),'Color',[0.9769 0.9839 0.0805]);

% Create semilogy
semilogy(X3,Y1,'DisplayName','$t(\dot{\epsilon}_\mathrm{min})$',...
    'LineWidth',2,...
    'Color',[1 0 0]);

% Create ylabel
ylabel('Stress [Pa]');

% Create xlabel
xlabel('$t-t_c$ [s]');

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes2,[-0.0160473127658922 0.8]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes2,[1e-5 1e4]);
hold(axes2,'off');
% Set the remaining axes properties
set(axes2,'YMinorTick','on','YScale','log');
% Create textarrow
annotation(figure1,'textarrow',[0.697666032303714,0.624543502659448],...
    [0.71830700061575,0.685541600877871],...
    'String',{'$t(\dot{\epsilon}_\mathrm{min})$'});

end














































