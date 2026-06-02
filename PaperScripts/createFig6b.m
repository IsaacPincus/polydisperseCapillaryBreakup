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
colours_MWDist = parula(Nm);
a0 = 6.3e-4;
X = 1;

opts = detectImportOptions("PS7MWDist.csv");
dataPS7 = readmatrix("PS7MWDist.csv", opts);
PS7MW = dataPS7(:,1);
PS7NormMass = dataPS7(:,2);

opts = detectImportOptions("PS16MWDist.csv");
dataPS16 = readmatrix("PS16MWDist.csv", opts);
PS16MW = dataPS16(:,1);
PS16NormMass = dataPS16(:,2);

% colours_relax_plots = [1,1,1;
%                        12,129,63;
%                        127,19,21;
%                        245,128,33]/255;

colours_relax_plots = flip(copper(4));

%% PS7_0 + PS16

opts = detectImportOptions(spreadsheetName, 'Sheet','PS16');
opts.SelectedVariableNames = [1:2,5:6,9:10,13:14,17:18,21:22]; 
% preview(spreadsheetName,opts)
M = readmatrix(spreadsheetName,opts);
t16 = M(2:end,1:2:11);
D16 = M(2:end,2:2:12);
% size16 = size(D16);

% project data onto same x-values
allPS_data = [PS7MW;PS16MW];
noPoints = 10000;
x = linspace(min(allPS_data), max(allPS_data),noPoints);
method = 'pchip';
yPS7 = interp1(PS7MW,PS7NormMass,x, method, 0);
yPS16 = interp1(PS16MW,PS16NormMass,x, method, 0);

% parameters for PS16
PS16_conc = [5,10,20,50,100,200];
PS7_conc = zeros(size(PS16_conc));
% PS16_conc = [0];
size16 = size(PS16_conc);

[~, t_c, tau_ep_min, tau_EC, ...
    aModel, epModel, tModel, stressModel, lambda_m] ...
        = get_outputs(PS16_conc, a0, eta_s, ...
            PS7_conc, yPS7, yPS16, vol_solvent, ...
            x, graphs_on, colours_MWDist, noPoints, ...
            Nm, nu, kB, T, NA, l, MW_monomer, cinf, ...
            Ns, K, alpha, chi, X, rho, cutoff, tmax, hsk, Ctau);

% %%

for ii = 1:size16(2)
% for ii = 1:1
    t = t16(:,ii);
    De = D16(:,ii);
    Dm = aModel{ii}/1e-3*2;
end

for ii = 1:size16(2)
% for ii = 1:1
    te = t16(:,ii);
    De = D16(:,ii);
    tavg = (te(2:end) + te(1:end-1))/2;
    Davg = (De(2:end) + De(1:end-1))/2;
    ep_exp = -2./Davg.*diff(De)./diff(te);
    % get relax times
    [ep_min_exp, index_ep_min] = min(ep_exp(tavg>0));
    tau_ep_min_exp(ii) = 2/(3*ep_min_exp);
    
    indices_EC = find((ep_exp<5*ep_min_exp)&(tavg>0));
    [xxx,yyy] = prepareCurveData(tavg(indices_EC),Davg(indices_EC));
    fitobject = fit(xxx, yyy, 'exp1',...
        'StartPoint', [0.1,-1/(3*tau_ep_min_exp(ii))]);
    tau_EC_fit = -1./(3*fitobject.b);
    tau_EC_exp(ii) = tau_EC_fit;

end

% %%
plotval = 1;
figure(59);
hold on
plot(PS16_conc, tau_EC, ':',...
    'Color',colours_relax_plots(plotval,:));
plot(PS16_conc, tau_EC_exp, 'o',...
    'Color',colours_relax_plots(plotval,:));
% plot(PS16_conc, tau_ep_min_exp, 's',...
%     'Color',colours_relax_plots(plotval,:));

axes1 = gca;
axes1.XScale = 'log';
axes1.YScale = 'log';

xlabel('$c_\mathrm{PS16}$ [ppm]')
ylabel('$\tau_\mathrm{EC}$ [s]')

figure(58);
hold on
plot(PS16_conc, tau_ep_min, '-',...
    'Color',colours_relax_plots(plotval,:));

axes1 = gca;
axes1.XScale = 'log';
axes1.YScale = 'log';

xlabel('$c_\mathrm{PS16}$ [ppm]')
ylabel('$\tau_\mathrm{EC}$ [s]')

tau_EC_ratio_PS16 = tau_EC/lambda_m(1);
tau_EC_exp_ratio_PS16 = tau_EC_exp/lambda_m(1);

%% PS7_200 + PS16

opts = detectImportOptions(spreadsheetName, 'Sheet','200ppmPS7+PS16');
opts.SelectedVariableNames = [1:2,5:6,9:10,13:14,17:18,21:22,25:26]; 
% preview(spreadsheetName,opts)
M = readmatrix(spreadsheetName,opts);
t16 = M(2:end,1:2:13);
D16 = M(2:end,2:2:14);
% size16 = size(D16);

% project data onto same x-values
allPS_data = [PS7MW;PS16MW];
noPoints = 10000;
x = linspace(min(allPS_data), max(allPS_data),noPoints);
method = 'pchip';
yPS7 = interp1(PS7MW,PS7NormMass,x, method, 0);
yPS16 = interp1(PS16MW,PS16NormMass,x, method, 0);

% parameters for PS16
PS16_conc = [0,5,10,20,50,100,200];
PS7_conc = 200*ones(size(PS16_conc));
% PS16_conc = [0];
size16 = size(PS16_conc);

[~, t_c, tau_ep_min, tau_EC, ...
    aModel, epModel, tModel, stressModel, lambda_m] ...
        = get_outputs(PS16_conc, a0, eta_s, ...
            PS7_conc, yPS7, yPS16, vol_solvent, ...
            x, graphs_on, colours_MWDist, noPoints, ...
            Nm, nu, kB, T, NA, l, MW_monomer, cinf, ...
            Ns, K, alpha, chi, X, rho, cutoff, tmax, hsk, Ctau);

for ii = 1:size16(2)
% for ii = 1:1
    t = t16(:,ii);
%     a = D16(:,ii)./D16(1,ii);
%     aM = aModel{ii}/D16(1,ii)/1e-3;
    De = D16(:,ii);
    Dm = aModel{ii}/1e-3*2;
end

for ii = 1:size16(2)
% for ii = 1:1
    te = t16(:,ii);
    De = D16(:,ii);
    tavg = (te(2:end) + te(1:end-1))/2;
    Davg = (De(2:end) + De(1:end-1))/2;
    ep_exp = -2./Davg.*diff(De)./diff(te);
    % get relax times
    [ep_min_exp, index_ep_min] = min(ep_exp(tavg>0&ep_exp>0));
    tau_ep_min_exp(ii) = 2/(3*ep_min_exp);
    
    indices_EC = find((ep_exp<5*ep_min_exp)&(tavg>0));
    [xxx,yyy] = prepareCurveData(tavg(indices_EC),Davg(indices_EC));
    fitobject = fit(xxx, yyy, 'exp1');
    tau_EC_fit = -1./(3*fitobject.b);
    tau_EC_exp(ii) = tau_EC_fit;

end

% %%
plotval = 2;
figure(59);
hold on
plot(PS16_conc, tau_EC, ':',...
    'Color',colours_relax_plots(plotval,:));
plot(PS16_conc, tau_EC_exp, 'o',...
    'Color',colours_relax_plots(plotval,:));
% plot(PS16_conc, tau_ep_min_exp, 's',...
%     'Color',colours_relax_plots(plotval,:));

figure(58);
hold on
plot(PS16_conc, tau_ep_min, '-',...
    'Color',colours_relax_plots(plotval,:));

tau_EC_ratio_PS7(1) = tau_EC(1)/lambda_m(2);
tau_EC_exp_ratio_PS7(1) = tau_EC_exp(1)/lambda_m(2);

%% PS7_500 + PS16

opts = detectImportOptions(spreadsheetName, 'Sheet','500ppmPS7+PS16');
opts.SelectedVariableNames = [1:2,5:6,9:10,13:14,17:18,21:22,25:26]; 
% preview(spreadsheetName,opts)
M = readmatrix(spreadsheetName,opts);
t16 = M(2:end,1:2:13);
D16 = M(2:end,2:2:14);
% size16 = size(D16);

% project data onto same x-values
allPS_data = [PS7MW;PS16MW];
noPoints = 10000;
x = linspace(min(allPS_data), max(allPS_data),noPoints);
method = 'pchip';
yPS7 = interp1(PS7MW,PS7NormMass,x, method, 0);
yPS16 = interp1(PS16MW,PS16NormMass,x, method, 0);

% parameters for PS16
PS16_conc = [0,5,10,20,50,100,200];
PS7_conc = 500*ones(size(PS16_conc));
% PS16_conc = [0];
size16 = size(PS16_conc);
X = 1;

[~, t_c, tau_ep_min, tau_EC, ...
    aModel, epModel, tModel, stressModel, lambda_m] ...
        = get_outputs(PS16_conc, a0, eta_s, ...
            PS7_conc, yPS7, yPS16, vol_solvent, ...
            x, graphs_on, colours_MWDist, noPoints, ...
            Nm, nu, kB, T, NA, l, MW_monomer, cinf, ...
            Ns, K, alpha, chi, X, rho, cutoff, tmax, hsk, Ctau);

for ii = 1:size16(2)
% for ii = 1:1
    t = t16(:,ii);
%     a = D16(:,ii)./D16(1,ii);
%     aM = aModel{ii}/D16(1,ii)/1e-3;
    De = D16(:,ii);
    Dm = aModel{ii}/1e-3*2;
end

for ii = 1:size16(2)
% for ii = 1:1
    te = t16(:,ii);
    De = D16(:,ii);
    tavg = (te(2:end) + te(1:end-1))/2;
    Davg = (De(2:end) + De(1:end-1))/2;
    ep_exp = -2./Davg.*diff(De)./diff(te);
    % get relax times
    [ep_min_exp, index_ep_min] = min(ep_exp(tavg>0));
    tau_ep_min_exp(ii) = 2/(3*ep_min_exp);
    
    indices_EC = find((ep_exp<5*ep_min_exp)&(tavg>0));
    [xxx,yyy] = prepareCurveData(tavg(indices_EC),Davg(indices_EC));
    fitobject = fit(xxx, yyy, 'exp1');
    tau_EC_fit = -1./(3*fitobject.b);
    tau_EC_exp(ii) = tau_EC_fit;
end

% %%
plotval = 3;
figure(59);
hold on
plot(PS16_conc, tau_EC, ':',...
    'Color',colours_relax_plots(plotval,:));
plot(PS16_conc, tau_EC_exp, 'o',...
    'Color',colours_relax_plots(plotval,:));
% plot(PS16_conc, tau_ep_min_exp, 's',...
%     'Color',colours_relax_plots(plotval,:));

figure(58);
hold on
plot(PS16_conc, tau_ep_min, '-',...
    'Color',colours_relax_plots(plotval,:));

tau_EC_ratio_PS7(2) = tau_EC(1)/lambda_m(2);
tau_EC_exp_ratio_PS7(2) = tau_EC_exp(1)/lambda_m(2);


%% PS7_1000 + PS16

opts = detectImportOptions(spreadsheetName, 'Sheet','1000ppmPS7+PS16');
opts.SelectedVariableNames = [1:2,5:6,9:10,13:14,17:18,21:22,25:26]; 
% preview(spreadsheetName,opts)
M = readmatrix(spreadsheetName,opts);
t16 = M(2:end,1:2:13);
D16 = M(2:end,2:2:14);
% size16 = size(D16);

% project data onto same x-values
allPS_data = [PS7MW;PS16MW];
noPoints = 10000;
x = linspace(min(allPS_data), max(allPS_data),noPoints);
method = 'pchip';
yPS7 = interp1(PS7MW,PS7NormMass,x, method, 0);
yPS16 = interp1(PS16MW,PS16NormMass,x, method, 0);

% parameters for PS16
PS16_conc = [0,5,10,20,50,100,200];
PS7_conc = 1000*ones(size(PS16_conc));
% PS16_conc = [0];
size16 = size(PS16_conc);
X = 1;

[~, t_c, tau_ep_min, tau_EC, ...
    aModel, epModel, tModel, stressModel, lambda_m] ...
        = get_outputs(PS16_conc, a0, eta_s, ...
            PS7_conc, yPS7, yPS16, vol_solvent, ...
            x, graphs_on, colours_MWDist, noPoints, ...
            Nm, nu, kB, T, NA, l, MW_monomer, cinf, ...
            Ns, K, alpha, chi, X, rho, cutoff, tmax, hsk, Ctau);

for ii = 1:size16(2)
% for ii = 1:1
    t = t16(:,ii);
%     a = D16(:,ii)./D16(1,ii);
%     aM = aModel{ii}/D16(1,ii)/1e-3;
    De = D16(:,ii);
    Dm = aModel{ii}/1e-3*2;
end

for ii = 1:size16(2)
% for ii = 1:1
    te = t16(:,ii);
    De = D16(:,ii);
    tavg = (te(2:end) + te(1:end-1))/2;
    Davg = (De(2:end) + De(1:end-1))/2;
    ep_exp = -2./Davg.*diff(De)./diff(te);
    % get relax times
    [ep_min_exp, index_ep_min] = min(ep_exp(tavg>0));
    tau_ep_min_exp(ii) = 2/(3*ep_min_exp);
    
    indices_EC = find((ep_exp<5*ep_min_exp)&(tavg>0));
    [xxx,yyy] = prepareCurveData(tavg(indices_EC),Davg(indices_EC));
    fitobject = fit(xxx, yyy, 'exp1');
    tau_EC_fit = -1./(3*fitobject.b);
    tau_EC_exp(ii) = tau_EC_fit;

end

% %%
plotval = 4;
figure(59);
hold on
plot(PS16_conc, tau_EC, ':',...
    'Color',colours_relax_plots(plotval,:));
plot(PS16_conc, tau_EC_exp, 'o',...
    'Color',colours_relax_plots(plotval,:));
% plot(PS16_conc, tau_ep_min_exp, 's',...
%     'Color',colours_relax_plots(plotval,:));


figure(58);
hold on
plot(PS16_conc, tau_ep_min, '-',...
    'Color',colours_relax_plots(plotval,:));
opts = detectImportOptions("Fig5_Hawardetal.csv");
data = readmatrix("Fig5_Hawardetal.csv", opts);
for ii = 1:4
    plot(data(:,ii*2-1),data(:,ii*2)/1000, 'p',...
        'Color',colours_relax_plots(ii,:), ...
        'MarkerFaceColor',colours_relax_plots(ii,:));
end

% combine them together into another plot
figure(60);
hold on
copyobj(figure(58).Children.Children, figure(60).Children);
copyobj(figure(59).Children.Children, figure(60).Children);
axes1 = gca;
axes1.XScale = 'log';
axes1.YScale = 'log';

xlabel('$c_\mathrm{PS16}$ [ppm]')
ylabel('$\tau_\mathrm{EC}$ [s]')

fig = figure(58);
fig.Position = [488,152.2,810,610];
legend({' ', ' ',' ', ' ', ...
    '$\qquad$ 0','$\qquad$ 200', '$\qquad$ 500','$\qquad$ 1000'},...
    'Position',[0.5379 0.2289 0.30303 0.19773],...
    'Orientation','vertical',...
    'NumColumns',2);
% Create textbox
annotation(fig,'textbox',...
    [0.6408 0.4269 0.0893 0.06356],...
    'String',{'Exp'});
% Create textbox
annotation(fig,'textbox',...
    [0.5375 0.4249 0.0864 0.063564],...
    'String',{'Sim'});
% Create textbox
annotation(fig,'textbox',...
    [0.74604 0.42495 0.15604 0.06356],...
    'String',{'$c_\mathrm{PS7}$'});
% Create textbox
annotation(fig,'textbox',...
    [0.54604 0.15495 0.5 0.06356],...
    'String',{'$\tau_\mathrm{EC}$ from $\min(\dot{\epsilon})$'});

fig = figure(59);
fig.Position = [488,152.2,810,610];
legend({' ', '$\qquad$ 0',...
    ' ', '$\qquad$ 200',...
    ' ', '$\qquad$ 500',...
    ' ', '$\qquad$ 1000'},...
    'Position',[0.5579 0.2289 0.30303 0.19773],...
    'Orientation','horizontal',...
    'NumColumns',2);
% Create textbox
annotation(fig,'textbox',...
    [0.6608 0.4269 0.0893 0.06356],...
    'String',{'Exp'});
% Create textbox
annotation(fig,'textbox',...
    [0.5575 0.4249 0.0864 0.063564],...
    'String',{'Sim'});
% Create textbox
annotation(fig,'textbox',...
    [0.76604 0.42495 0.15604 0.06356],...
    'String',{'$c_\mathrm{PS7}$'});
% Create textbox
annotation(fig,'textbox',...
    [0.56604 0.15495 0.5 0.06356],...
    'String',{'$\tau_\mathrm{EC}$ from $D \sim e^{-t/3\tau_\mathrm{EC}}$'});

% fig = figure(60);
% fig.Position = [488,152.2,810,610];
% legend({'$\tau_\mathrm{EC}$ from $\min(\dot{\epsilon})$, sim','','','',...
%     '$\tau_\mathrm{EC}$ from $\min(\dot{\epsilon})$, exp','','','',...
%     '$\tau_\mathrm{EC}$ from $D \sim e^{-t/3\tau_\mathrm{EC}}$, sim',...
%     '$\tau_\mathrm{EC}$ from $D \sim e^{-t/3\tau_\mathrm{EC}}$, exp',...
%     '','','','','',''})

tau_EC_ratio_PS7(3) = tau_EC(1)/lambda_m(2);
tau_EC_exp_ratio_PS7(3) = tau_EC_exp(1)/lambda_m(2);

%% tau/tau_EC vs c/c*

cstar_PS7 = 3.92E+03; % wppm
cstar_PS16 = 2.59E+03; % wppm

cPS7 = [200, 500, 1000];
cPS16 = [5,10,20,50,100,200];

cr7 = cPS7/cstar_PS7;
cr16 = cPS16/cstar_PS16;

figure();
hold on
plot(cr16, tau_EC_ratio_PS16, 'rs', 'MarkerFaceColor','r', ...
    'DisplayName','PS16 sim', ...
    'MarkerSize', 12, 'LineWidth', 2);
plot(cr16, tau_EC_exp_ratio_PS16, 'rs', 'MarkerFaceColor','none', ...
    'DisplayName','PS16 exp', ...
    'MarkerSize', 12, 'LineWidth', 2);
plot(cr7, tau_EC_ratio_PS7, 'bo', 'MarkerFaceColor','b', ...
    'DisplayName','PS7 sim', ...
    'MarkerSize', 12, 'LineWidth', 2);
plot(cr7, tau_EC_exp_ratio_PS7, 'bo', 'MarkerFaceColor','none', ...
    'DisplayName','PS7 exp', ...
    'MarkerSize', 12, 'LineWidth', 2);

leg = legend;
leg.NumColumns = 2;

xlabel('$c/c^*$')
ylabel('$\tau_\mathrm{EC}/\tau_Z$')

axes1 = gca;
axes1.XScale = 'log';

xlim([1e-3, 4e-1]);

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
annotation(fig,'line',[0.481719367588934 0.454051383399211],...
    [0.297165137614677 0.224770642201833],'Color',[0 0 1],'LineWidth',0.5,...
    'LineStyle','--');

% Create line
annotation(fig,'line',[0.350790513833992 0.323122529644269],...
    [0.377768020969855 0.305373525557011],'Color',[0 0 1],'LineWidth',0.5,...
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
        
        [a,~,~,~,stress,ep,~,t,~,~,~] = ...
            filamentThinningFENEPM_inertia_CDD(chi, X, mu, rho, ...
            Nm, Ns, L, Nk, hsk, lambda_ms, g, a0, cutoff, tmax, phieq);

        % [a,~,~,~,stress,ep,~,t,~,~,~] = ...
        %     filamentThinningFENEPM_inertia(chi, X, mu, rho, ...
        %     Nm, Ns, L, lambda_ms, g, a0, cutoff, tmax);
        
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