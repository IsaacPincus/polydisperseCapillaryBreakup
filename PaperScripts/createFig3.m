clear variables;
close all;
warning('off','MATLAB:nearlySingularMatrix');

datapath = "../Data";
addpath(datapath)
addpath("../Functions")

spreadsheetName = 'ThinningCurves_PS7_PS16_andMixtures.xlsx';

% PS constants
graphs_on = false;
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
Nm = 1;
Ns = 30;
hstar = 0.267;
hsk = 0.027;
phieq = 0;
Ctau = 1;
a0 = 6.3e-4;
X = 1;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now do subplots
% clear variables

PS7_conc = 1000;
PS16_conc = 500;

colours_MWDist = parula(3);
colours_MWDist = colours_MWDist([1,2],:);
graphs_on = false;

opts = detectImportOptions("PS7MWDist.csv");
dataPS7 = readmatrix("PS7MWDist.csv", opts);
PS7MW = dataPS7(:,1);
PS7NormMass = dataPS7(:,2);

opts = detectImportOptions("PS16MWDist.csv");
dataPS16 = readmatrix("PS16MWDist.csv", opts);
PS16MW = dataPS16(:,1);
PS16NormMass = dataPS16(:,2);


%% fit individual, two modes

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

Nm = 2;
mol_polymer_density = mol_polymer_densityPS7 + mol_polymer_densityPS16;
mol_frac_PS7 = mol_polymer_densityPS7/mol_polymer_density;
n = [nPS16*(1-mol_frac_PS7), nPS7*mol_frac_PS7];
M = [MPS16, MPS7];

% get densities etc
density = mol_polymer_density; % mol/m^3, total mol density of polymer
g = n'*density*kB*T*NA;

Ltrue = 2*l*M/MW_monomer;
R0 = sqrt(2*l^2.*M/MW_monomer*cinf);
Rg = R0/sqrt(6);
Nk = Ltrue.^2./(2*cinf*l^2.*M/MW_monomer);
Nk = Nk';
% L = Ltrue'./R0';
L = sqrt(3*((Ltrue'./R0').^2-1));

jvals = 1:Ns;
%     hstar = 0;
sigma = -1.40*hstar^0.78;
b = 1-1.66*hstar^0.78;
intrinsic_visc = K*M.^alpha;
lambda_m = M.*intrinsic_visc*eta_s/(NA*kB*T);
%     lambda_m = M.^1.5.*intrinsic_visc*eta_s/(NA*kB*T);
lambda_ms = lambda_m'./jvals.^(2+sigma);

%% solve for stresses

% [a,~,~,~,stress,ep,~,t,~,~,~] = ...
%     filamentThinningFENEPM_inertia_CDD(chi, X, mu, rho, ...
%     Nm, Ns, L, Nk, hsk, lambda_ms, g, a0, cutoff, tmax, phieq);

[a,Az,Ar,f,stress,ep,Wi,t,te,ye,ie] = ...
    filamentThinningFENEPM_inertia(chi, X, mu, rho, ...
    Nm, Ns, L, lambda_ms, g, a0, cutoff, tmax);

local_ep_maxes = find(islocalmax(ep)==1);
t_max_ep_all = t(islocalmax(ep)==1);
t_c = t_max_ep_all(1);

ep_min_array = islocalmin(ep)==1;
ep_min_all = ep(ep_min_array&(t>t_c));
ep_min = min(ep_min_all);
t_ep_min = min(t(ep_min_array&(t>t_c)));

tau_ep_min = 2/(3*ep_min);

% indices_EC = find((ep<1.2*ep_min)&(t>t_c));
indices_EC = find((ep<1.2*ep_min)&(t>t_c));
[xxx,yyy] = prepareCurveData(t(indices_EC)-t_c,a(indices_EC));
fitobject = fit(xxx, yyy, 'exp1');
tau_EC_fit = -1./(3*fitobject.b);
tau_EC = tau_EC_fit;

%     t_c(ii) = 0;

aModel = a;
epModel = ep;
tModel = t;
stressModel = stress;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting fits
% fig = figure(59);
% axes1 = axes(fig, 'Position',...
    % [0.28365612 0.20061 0.27958 0.28955]);
fig = figure();
axes1 = gca;
hold on
% axes1.FontSize = 13;
axes1.YScale = 'log';
tm = t-t_c;
fitPlot = fitobject;
fitPlot.a = fitPlot.a*1e3;
plot(fitPlot)
plot(tm, a*1e3, 'k');
plot(tm(indices_EC), a(indices_EC)*1e3, 'or');
legend({'Exponential Fit', 'Simulation', 'Points used for fit'})
xlabel('$t-t_c$ [s]');
ylabel('$a$ [mm]');
xlim([-0.015,1.2]);
ylim([1e-6,3e-4]*1e3);

% Create line
annotation(fig,'line',[0.282219085262563,0.25455110107284],...
    [0.722529925731762,0.650135430318918],'Color',[0 0 1],'LineWidth',0.5,...
    'LineStyle','--');

% Create line
annotation(fig,'line',[0.577433653303219,0.549765669113496],...
    [0.491450851900393,0.419056356487549],'Color',[0 0 1],'LineWidth',0.5,...
    'LineStyle','--');

% fig = figure(58);
% axes1 = axes(fig, 'Position',...
%     [0.257727272727274 0.208514637836674 0.267727272727273 0.294761902137113]);
fig = figure();
axes1 = gca;
hold on
axes1.YScale = 'log';
% axes1.FontSize = 13;
tm = t-t_c;
plot(tm, ep, 'k');
plot(t_ep_min-t_c, ep_min, 'or');
plot(tm, ones(size(ep))*ep_min, 'b--', 'linewidth', 0.5);
plot(tm, ones(size(ep))*ep_min*1.2, 'b--', 'linewidth', 0.5);
legend({'Simulation', 'Min $\dot{\epsilon}$'})
xlabel('$t-t_c$ [s]');
ylabel('$\dot{\epsilon}$ [s$^{-1}$]');
xlim([-0.015,1.2]);
ylim([4,1000]);



%% functions

function [a0, t_c, tau_ep_min, tau_EC,...
    aModel, epModel, tModel, stressModel, lambda_m] ...
    = get_outputs(PS16_conc, a0, eta_s, ...
    PS7_conc, yPS7, yPS16, vol_solvent, ...
    x, graphs_on, colours_MWDist, noPoints, ...
    Nm, nu, kB, T, NA, l, MW_monomer, cinf, ...
    Ns, K, alpha, chi, X, rho, cutoff, tmax, hsk, Ctau)

    for ii = 1:length(PS16_conc)
        % for ii = 1:1
        % S = 1;
        S = 2.36;
        
%         a0 = 3.5e-4;
        mu = eta_s;
        %     texp = t16(:,ii);
        %     Dexp = D16(:,ii)*1e-3;

%         graphs_on = true;

        Nm = 1;
        % fit in two stages
        mass_polymer = PS7_conc(ii);
        weight_fraction_PS16 = 0;
        y = yPS7*(1-weight_fraction_PS16)+yPS16*weight_fraction_PS16;
        density_mass_PS = mass_polymer/vol_solvent;
        maxMW = max(x);
        [nPS7, MPS7, mol_polymer_densityPS7] = ...
            splitMWDist(graphs_on, colours_MWDist,...
            x, y, noPoints, maxMW, Nm, nu, density_mass_PS);
        
        mass_polymer = PS16_conc(ii);
        weight_fraction_PS16 = 1;
        y = yPS7*(1-weight_fraction_PS16)+yPS16*weight_fraction_PS16;
        density_mass_PS = mass_polymer/vol_solvent;
        maxMW = max(x);
        [nPS16, MPS16, mol_polymer_densityPS16] = ...
            splitMWDist(graphs_on, colours_MWDist,...
            x, y, noPoints, maxMW, Nm, nu, density_mass_PS);

        mol_polymer_density = mol_polymer_densityPS7 + mol_polymer_densityPS16;
        mol_frac_PS7 = mol_polymer_densityPS7/mol_polymer_density;
        n = [nPS16*(1-mol_frac_PS7), nPS7*mol_frac_PS7];
        M = [MPS16, MPS7];
        Nm = 2;
        
        % get densities etc
        density = mol_polymer_density; % mol/m^3, total mol density of polymer
        g = n'*density*kB*T*NA;
        % change g by some function, see what happens
%         g = g.*(1+5e-8*M');
        
        Ltrue = 2*l*M/MW_monomer;
        R0 = sqrt(2*l^2.*M/MW_monomer*cinf);
        Rg = R0/sqrt(6);
        Nk = Ltrue.^2./(2*cinf*l^2.*M/MW_monomer);
        Nk = Nk';
        % L = Ltrue'./R0';
        L = sqrt(3*((Ltrue'./R0').^2-1));
        % L = sqrt(3*((Ltrue'./R0').^2-1));

        phieq = 0;
        
        jvals = 1:Ns;
        hstar = 0.2671;
        % hstar = 0.029;
        %     hstar = 0;
        sigma = -1.40*hstar^0.78;
        b = 1-1.66*hstar^0.78;
        intrinsic_visc = K*M.^alpha;
        lambda_m = 1/S*M.*intrinsic_visc*eta_s/(NA*kB*T)*Ctau;
        %     lambda_m = M.^1.5.*intrinsic_visc*eta_s/(NA*kB*T);
        lambda_ms = lambda_m'./jvals.^(2+sigma);
        
        % [a,~,~,~,stress,ep,~,t,~,~,~] = ...
        %     filamentThinningFENEPM_inertia_CDD(chi, X, mu, rho, ...
        %     Nm, Ns, L, Nk, hsk, lambda_ms, g, a0, cutoff, tmax, phieq);

        [a,~,~,~,stress,ep,~,t,~,~,~] = ...
            filamentThinningFENEPM_inertia(chi, X, mu, rho, ...
            Nm, Ns, L, lambda_ms, g, a0, cutoff, tmax);
        
        local_ep_maxes = find(islocalmax(ep)==1);
        t_max_ep_all = t(islocalmax(ep)==1);
        t_c(ii) = t_max_ep_all(1);
        
        ep_min_array = islocalmin(ep)==1;
        ep_min_all = ep(ep_min_array&(t>t_c(ii)));
        ep_min = min(ep_min_all);
        
        tau_ep_min(ii) = 2/(3*ep_min);
        
        indices_EC = find((ep<1.2*ep_min)&(t>t_c(ii)));
        [xxx,yyy] = prepareCurveData(t(indices_EC),a(indices_EC));
        fitobject = fit(xxx, yyy, 'exp1');
        tau_EC_fit = -1./(3*fitobject.b);
        tau_EC(ii) = tau_EC_fit;
        
        %     t_c(ii) = 0;
        
        aModel{ii} = a;
        epModel{ii} = ep;
        tModel{ii} = t;
        stressModel{ii} = stress;
    end
end