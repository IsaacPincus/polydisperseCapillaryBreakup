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
Ns = 30;
hstar = 0.267;
hsk = 0.027;
% hsk = 0.0;
phieq = 0;
% Ctau = 1;
Uetatau = 2.36;
a0 = 6.3e-4;
X = 1;
colours_MWDist = parula(6);
colours_MWDist = colours_MWDist([1,5],:);
graphs_on = false;
PS7_conc = 1000;
PS16_conc = 100;

opts = detectImportOptions("PS7MWDist.csv");
dataPS7 = readmatrix("PS7MWDist.csv", opts);
PS7MW = dataPS7(:,1);
PS7NormMass = dataPS7(:,2);

opts = detectImportOptions("PS16MWDist.csv");
dataPS16 = readmatrix("PS16MWDist.csv", opts);
PS16MW = dataPS16(:,1);
PS16NormMass = dataPS16(:,2);


% %% fit individual, two modes

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
lambda_m = 1/Uetatau*M.*intrinsic_visc*eta_s/(NA*kB*T);
%     lambda_m = M.^1.5.*intrinsic_visc*eta_s/(NA*kB*T);
lambda_ms = lambda_m'./jvals.^(2+sigma);

% %% solve for stresses

[a,~,~,~,stress,ep,Wi,t,~,~,~] = ...
    filamentThinningFENEPM_inertia_CDD(chi, X, mu, rho, ...
    Nm, Ns, L, Nk, hsk, lambda_ms, g, a0, cutoff, tmax, phieq);

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

% %% plotting individual for last run

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% contributions to stress
fig = figure();
fig.Position = [488,152.2,810,610];
hold on

% capillary stress
plot(tModel-t_c, chi./aModel, 'k', 'linewidth', 2);

xlabel('$t-t_c$ [s]')
% ylabel('$g_i f_i (A_{z,i} - A_{r,i})$')
ylabel('stress [Pa]')

axes1 = gca;
axes1.YScale = 'log';
% axes1.XScale = 'log';

% polymer stresses
for pp = 1:Nm
    plot(tModel-t_c, stress(:,pp), '-', 'Color',colours_MWDist(pp,:));
    % heightIndex = find(tModel-t_c>2*lambda_ms(pp),1);
    % plot(2*[lambda_ms(pp),lambda_ms(pp)], ...
    %       [1e-5,stress(heightIndex,pp)], ':',...
    %       'HandleVisibility','off',...
    %       'Color',colours_MWDist(pp,:))
end

% intertial and viscous stresses
plot(tModel-t_c, 3*eta_s*ep + 1/8*rho*ep.^2.*aModel.^2, ...
    'k-.', 'linewidth', 1.5);

l1 = legend({'Capillary Stress',...
    sprintf('PS16 Polymer Stress'),...
    sprintf('PS7 Polymer Stress'),...
    'Viscous and Intertial Stresses'});

xlim([tModel(1)-t_c,tModel(end)-t_c])
ylim([8e-2,1e5])

copyobj(fig.Children(2), fig);
axes2 = fig.Children(3);
axes2.Position = [0.465454545454546,0.441181102362205,0.218,0.23];
axes2.XTick = [];
axes2.YTick = [];
xlabel('')
ylabel('')
xlim([tModel(1)-t_c,0.02])
ylim([3,700])
fig.Children = flip(fig.Children);

l1.Position = [0.162,0.719,0.388,0.174];

% Create arrow
annotation(fig,'arrow',[0.464426877470356,0.189723320158103],...
    [0.606882281916595,0.537419634472297]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extension rate
fig = figure();
fig.Position = [488,152.2,810,610];
hold on

names = {'PS16', 'PS7'};

for pp = 1:Nm
    plot(tModel-t_c, Wi(pp,:),...
        '-', 'Color',colours_MWDist(pp,:), ...
        'DisplayName',names{pp});
end

xlabel('$t-t_c$ [s]')
ylabel('$W\!i_\mathrm{eff}$')
legend('Location','northeast')

axes1 = gca;
axes1.YScale = 'log';
% axes1.XScale = 'log';

xlim([tModel(1)-t_c,tModel(end)-t_c])
ylim([5e-2,1e2])

plot(tModel-t_c, ones(size(tModel))*0.5, 'k--', ...
    'DisplayName','$W\!i_\mathrm{eff} = 0.5$ (coil-stretch)')

copyobj(fig.Children(2), fig);
axes2 = fig.Children(1);
axes2.Position = [0.289555555555555,0.654241469816274,0.218,0.23];
% axes2.XTick = [];
% axes2.YTick = [];
xlabel(axes2, '')
ylabel(axes2, '')
xlim(axes2, [tModel(1)-t_c,0.01])
ylim(axes2, [2e-1,9e1])

% fig.Children = flip(fig.Children);
% Create arrow
annotation(fig,'arrow',[0.280740740740741,0.184782608695652],...
    [0.801837270341207,0.742782152230971]);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables;

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
Nm = 30;
Ns = 1;
hstar = 0.267;
hsk = 0.027;
% hsk = 0.0;
phieq = 0;
% Ctau = 1;
Uetatau = 2.36;
a0 = 6.3e-4;
PS7_conc = 1000;
PS16_conc = 100;
mu = eta_s;
X = 1;
colours_MWDist = parula(Nm);
graphs_on = false;

opts = detectImportOptions("PS7MWDist.csv");
dataPS7 = readmatrix("PS7MWDist.csv", opts);
PS7MW = dataPS7(:,1);
PS7NormMass = dataPS7(:,2);

opts = detectImportOptions("PS16MWDist.csv");
dataPS16 = readmatrix("PS16MWDist.csv", opts);
PS16MW = dataPS16(:,1);
PS16NormMass = dataPS16(:,2);

%% fit together, multiple modes

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

%% solve for stresses
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
lambda_m = 1/Uetatau*M.*intrinsic_visc*eta_s/(NA*kB*T);
%     lambda_m = M.^1.5.*intrinsic_visc*eta_s/(NA*kB*T);
lambda_ms = lambda_m'./jvals.^(2+sigma);

[a,~,~,~,stress,ep,Wi,t,~,~,~] = ...
    filamentThinningFENEPM_inertia_CDD(chi, X, mu, rho, ...
    Nm, Ns, L, Nk, hsk, lambda_ms, g, a0, cutoff, tmax, phieq);

local_ep_maxes = find(islocalmax(ep)==1);
t_max_ep_all = t(islocalmax(ep)==1);
t_c = t_max_ep_all(1);

ep_min_array = islocalmin(ep)==1;
ep_min_all = ep(ep_min_array&(t>t_c));
ep_min = min(ep_min_all);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% contributions to stress
fig = figure();
fig.Position = [488,152.2,810,610];
hold on

% capillary stress
plot(tModel-t_c, chi./aModel, 'k', 'linewidth', 2);

% intertial and viscous stresses
plot(tModel-t_c, 3*eta_s*ep + 1/8*rho*ep.^2.*aModel.^2, ...
    'k-.', 'linewidth', 1.5);

xlabel('$t-t_c$ [s]')
ylabel('Stress [Pa]')

axes1 = gca;
axes1.YScale = 'log';
% axes1.XScale = 'log';

% polymer stresses
for pp = 1:Nm
    p1 = plot(tModel-t_c, stress(:,pp), '-', 'Color',colours_MWDist(pp,:));
%     heightIndex = find(tModel-t_c>2*lambda_ms(pp),1);
%     plot(2*[lambda_ms(pp),lambda_ms(pp)], ...
%           [1e-5,stress(heightIndex,pp)], ':',...
%           'HandleVisibility','off',...
%           'Color',colours_MWDist(pp,:))
end

l1 = legend({'Capillary Stress',...
    'Viscous and Intertial Stresses'});

xlim([tModel(1)-t_c,tModel(end)-t_c])
ylim([5e-5,1e5])

lines = fig.Children(2).Children;

lines(16).LineStyle = '--';
lines(16).LineWidth = 5;
uistack(lines(16),'top');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extension rate
fig = figure();
fig.Position = [488,152.2,810,610];
hold on

for pp = 1:Nm
    plot(tModel-t_c, Wi(pp,:),...
        '-', 'Color',colours_MWDist(pp,:));
end

xlabel('$t-t_c$ [s]')
ylabel('$W\!i_{\mathrm{eff}, i}$')
% legend('Location','northeast')

axes1 = gca;
axes1.YScale = 'log';
% axes1.XScale = 'log';

xlim([tModel(1)-t_c,tModel(end)-t_c])
ylim([1e-2,1e2])

plot(tModel-t_c, ones(size(tModel))*0.5, 'k--', ...
    'DisplayName','$W\!i_\mathrm{eff} = 0.5$ (coil-stretch)')

lines = fig.Children(1).Children;

lines(16).LineStyle = '--';
lines(16).LineWidth = 5;
uistack(lines(16),'top');

copyobj(fig.Children(1), fig);
axes2 = fig.Children(1);
axes2.Position = [0.342888888888889,0.641118110236221,0.218,0.23];
% axes2.XTick = [];
% axes2.YTick = [];
xlabel(axes2, '')
ylabel(axes2, '')
xlim(axes2, [tModel(1)-t_c,0.01])
ylim(axes2, [5e-2,9e1])

% fig.Children = flip(fig.Children);
% Create arrow
annotation(fig,'arrow',[0.33498023715415,0.184782608695652],...
    [0.784776902887139,0.742782152230971]);







































