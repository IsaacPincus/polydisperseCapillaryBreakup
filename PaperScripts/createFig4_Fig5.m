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
% eta_s = 100; % Pa.s, solvent viscosity
NA = 6.02e23;  % 1/mol, Avagadro's constant
kB = 1.38e-23;  % J/k, Boltzmann constant
T = 298;        % K, temperature
l = 1.54e-10;   % m, PS bond length
cinf = 9.7;     % PS characteristic ratio
MW_monomer = 104; % g/mol, PS monomer weight
nu = 0.5;
chi = 3.1e-2; % N/m, DOP maybe 0.015 N/m?
rho = 985; % kg/m^3, fluid density
alpha_C1 = 0.25;
C1 = 9/(4*alpha_C1^3);
rho = rho*C1*8;
tmax = 100; % s 
cutoff = 1e-3;
Nm = 1;
Ns = 10;
hsk = 0.027;
Ctau = 1;
Uetatau = 2.36;
colours_MWDist = parula(Nm);
% a0 = 6.2e-4;
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
M{1} = readmatrix(spreadsheetName,opts);
t16 = M{1}(2:end,1:2:11);
D16 = M{1}(2:end,2:2:12);
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

figD = figure();
figD.Position = [488,152.2,810,610];
hold on
colours = winter(size16(2)+1);
for ii = 1:size16(2)
% for ii = 1:1
    t = t16(:,ii);
    ae = D16(:,ii)/2;
    am = aModel{ii}/1e-3;
    plot(t,ae, 's:', 'Color', colours(ii,:),...
        'DisplayName',sprintf('$c_\\mathrm{PS16} = %d$ [wppm]', PS16_conc(ii)))
    plot(tModel{ii}-t_c(ii), am, '-', 'linewidth', 2, ...
        'Color',colours(ii,:), 'HandleVisibility','off');
end
axes1 = gca;
axes1.YScale = 'log';
xlabel('$t-t_c$ [s]');
% ylabel('$D/D_0$');
ylabel('$a$ [mm]');
l1 = legend;
l1.Position = [0.54485,0.5191,0.3131,0.2941];

xlim(axes1,[-0.05 0.65]);
ylim(axes1,[5e-3 1.005]/2);

% % Create textbox
% annotation('textbox', ...
%     [0.38399 0.56304 0.273818 0.353342],...
%     'String',...
%     ["$\tau_i = \frac{M_i [\eta] \eta_s }{N_A k_\mathrm{B} T}$",...
%     "$[\eta] = K M^\alpha$",...
%     sprintf("$\\eta_s = %0.4g$ [Pa.s]", eta_s),...
%     sprintf("$K = %0.4g$", K),...
%     sprintf("$\\alpha = %0.4g$", alpha),...
%     sprintf("$\\chi = %0.4g$ [N m$^{-1}$]", chi),...
%     sprintf("$D_0 = %0.4g$ [mm]", a0*1e3*2),...
%     ])

% Create textbox
annotation('textbox',...
    [0.586,0.825,0.251,0.070],...
    'String',{'$c_\mathrm{PS7} = 0$ [wppm]'},...
    'FitBoxToText','off');

% %%
figep = figure();
figep.Position = [488,152.2,810,610];
hold on
colours = winter(size16(2)+1);
for ii = 1:size16(2)
% for ii = 1:1
    te = t16(:,ii);
    ae = D16(:,ii)/2;
    tavg = (te(2:end) + te(1:end-1))/2;
    aavg = (ae(2:end) + ae(1:end-1))/2;
    ep_exp = -2./aavg.*diff(ae)./diff(te);
    % get relax times
    [ep_min_exp, index_ep_min] = min(ep_exp(tavg>0));
    tau_ep_min_exp(ii) = 2/(3*ep_min_exp);

    indices_EC = find((ep_exp<5*ep_min_exp)&(tavg>0));
    [xxx,yyy] = prepareCurveData(tavg(indices_EC),aavg(indices_EC));
    fitobject = fit(xxx, yyy, 'exp1',...
        'StartPoint', [0.1,-1/(3*tau_ep_min_exp(ii))]);
    tau_EC_fit = -1./(3*fitobject.b);
    tau_EC_exp(ii) = tau_EC_fit;

    plot(tavg,ep_exp, 's', 'Color', colours(ii,:),...
        'DisplayName',sprintf('%d ppm PS16', PS16_conc(ii)))
    plot(tModel{ii}-t_c(ii), epModel{ii}, '-', 'linewidth', 2, ...
        'Color',colours(ii,:), 'HandleVisibility','off');
end
axes1 = gca;
axes1.YScale = 'log';
xlabel('$t-t_c$ [s]');
% ylabel('$D/D_0$');
ylabel('$\dot{\epsilon}$ [1/s]');
% legend
xlim(axes1,[-0.05 0.65]);
ylim(axes1,[4,2000]);

tau_EC_ratio_PS16 = tau_EC/lambda_m(1);
tau_EC_exp_ratio_PS16 = tau_EC_exp/lambda_m(1);

M_implied = (tau_EC_exp/(K*eta_s)*(NA*kB*T*Uetatau)).^(2/3);

%% PS16_0 + PS7

clear tau_EC tau_EC_exp

opts = detectImportOptions(spreadsheetName, 'Sheet','200ppmPS7+PS16');
opts.SelectedVariableNames = [1:2]; 
M{1} = readmatrix(spreadsheetName,opts);
opts = detectImportOptions(spreadsheetName, 'Sheet','500ppmPS7+PS16');
opts.SelectedVariableNames = [1:2]; 
M{2} = readmatrix(spreadsheetName,opts);
opts = detectImportOptions(spreadsheetName, 'Sheet','1000ppmPS7+PS16');
opts.SelectedVariableNames = [1:2]; 
M{3} = readmatrix(spreadsheetName,opts);

% project data onto same x-values
allPS_data = [PS7MW;PS16MW];
noPoints = 10000;
x = linspace(min(allPS_data), max(allPS_data),noPoints);
method = 'pchip';
yPS7 = interp1(PS7MW,PS7NormMass,x, method, 0);
yPS16 = interp1(PS16MW,PS16NormMass,x, method, 0);

% parameters for PS16
PS16_conc = [0,0,0];
PS7_conc = [200, 500, 1000];
% PS16_conc = [0];
size7 = size(PS7_conc);
X = 1;

[~, t_c, tau_ep_min, tau_EC, ...
    aModel, epModel, tModel, stressModel, lambda_m] ...
        = get_outputs(PS16_conc, a0, eta_s, ...
            PS7_conc, yPS7, yPS16, vol_solvent, ...
            x, graphs_on, colours_MWDist, noPoints, ...
            Nm, nu, kB, T, NA, l, MW_monomer, cinf, ...
            Ns, K, alpha, chi, X, rho, cutoff, tmax, hsk, Ctau);

% %%
figD = figure();
figD.Position = [488,152.2,810,610];
hold on
colours = copper(size7(2)+1);
colours = flipud(colours(1:end-1,:));
for ii = 1:size7(2)
% for ii = 1:1
    t = M{ii}(:,1);
%     a = D16(:,ii)./D16(1,ii);
%     aM = aModel{ii}/D16(1,ii)/1e-3;
    ae = M{ii}(:,2)/2;
    am = aModel{ii}/1e-3;
    plot(t,ae, 's:', 'Color', colours(ii,:),...
        'DisplayName',sprintf('$c_\\mathrm{PS7} = %d$ [wppm]', PS7_conc(ii)))
    plot(tModel{ii}-t_c(ii), am, '-', 'linewidth', 2, ...
        'Color',colours(ii,:), 'HandleVisibility','off');
end
axes1 = gca;
axes1.YScale = 'log';
xlabel('$t-t_c$ [s]');
% ylabel('$D/D_0$');
ylabel('$a$ [mm]');
l1 = legend;
l1.Position = [0.551,0.627,0.32,0.19];

xlim(axes1,[-0.05 0.65]);
ylim(axes1,[5e-3 1.005]/2);

% Create textbox
annotation('textbox',...
    [0.586,0.824,0.30,0.07],...
    'String',{'$c_\mathrm{PS16} = 0$ [wppm]'},...
    'FitBoxToText','off');

% %%
figep = figure();
figep.Position = [488,152.2,810,610];
hold on
% colours = flipud(copper(size7(2)+1));
for ii = 1:size7(2)
% for ii = 1:1
    te = M{ii}(:,1);
    ae = M{ii}(:,2)/2;
    tavg = (te(2:end) + te(1:end-1))/2;
    aavg = (ae(2:end) + ae(1:end-1))/2;
    ep_exp = -2./aavg.*diff(ae)./diff(te);
    % get relax times
    [ep_min_exp, index_ep_min] = min(ep_exp(tavg>0));
    tau_ep_min_exp(ii) = 2/(3*ep_min_exp);

    indices_EC = find((ep_exp<3*ep_min_exp)&(tavg>0));
    [xxx,yyy] = prepareCurveData(tavg(indices_EC),aavg(indices_EC));
    fitobject = fit(xxx, yyy, 'exp1');
    tau_EC_fit = -1./(3*fitobject.b);
    tau_EC_exp(ii) = tau_EC_fit;

    plot(tavg,ep_exp, 's', 'Color', colours(ii,:),...
        'DisplayName',sprintf('$c_\\mathrm{PS16} = %d$ [wppm]', PS16_conc(ii)))
    plot(tModel{ii}-t_c(ii), epModel{ii}, '-', 'linewidth', 2, ...
        'Color',colours(ii,:), 'HandleVisibility','off');
end
axes1 = gca;
axes1.YScale = 'log';
xlabel('$t-t_c$ [s]');
% ylabel('$D/D_0$');
ylabel('$\dot{\epsilon}$ [1/s]');

xlim(axes1,[-0.05 0.65]);
ylim(axes1,[4 2000]);

tau_EC_ratio_PS7 = tau_EC/lambda_m(2);
tau_EC_exp_ratio_PS7 = tau_EC_exp/lambda_m(2);

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

        phieq = 0;
        Uetatau = 2.36;
        S = Uetatau;
        
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
%         L = Ltrue'./R0';
        L = sqrt(3*((Ltrue'./R0').^2-1));
        
        jvals = 1:Ns;
        hstar = 0.2671;
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