clear variables;
close all;
warning('off','MATLAB:nearlySingularMatrix');

datapath = "../Data";
addpath(datapath)
addpath("../Functions")

% PS constants
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
% rho = rho*C1*8;
tmax = 100; % s
cutoff = 1e-4;
Ns = 50;
% hstar = 0.2671;
hstar = 0.26;
% hstar = 0.1;
mu = eta_s;
a0 = 6.3e-4;  % initial filament radius, m
X = 1;
Nm = 1;
S = 2.36;
colours_MWDist = parula(Nm+1);
colours_MWDist = colours_MWDist(1:Nm,:);
graphs_on = false;


opts = detectImportOptions("PS7MWDist.csv");
dataPS7 = readmatrix("PS7MWDist.csv", opts);
PS7MW = dataPS7(:,1);
PS7NormMass = dataPS7(:,2);

opts = detectImportOptions("PS16MWDist.csv");
dataPS16 = readmatrix("PS16MWDist.csv", opts);
PS16MW = dataPS16(:,1);
PS16NormMass = dataPS16(:,2);


%% 1000 ppm PS7, no PS16

PS7_conc = 1000;
PS16_conc = 0;

% project data onto same x-values
allPS_data = [PS7MW;PS16MW];
noPoints = 10000;
x = linspace(min(allPS_data), max(allPS_data),noPoints);
method = 'pchip';
yPS7 = interp1(PS7MW,PS7NormMass,x, method, 0);
yPS16 = interp1(PS16MW,PS16NormMass,x, method, 0);

mass_polymer = PS7_conc;

weight_fraction_PS16 = 0/mass_polymer;
y = yPS7*(1-weight_fraction_PS16)+yPS16*weight_fraction_PS16;
density_mass_PS = mass_polymer/vol_solvent;
maxMW = max(x);
[n, M, mol_polymer_density] = ...
    splitMWDist(graphs_on, colours_MWDist,...
    x, y, noPoints, maxMW, Nm, nu, density_mass_PS);

%%
% get densities etc
density = mol_polymer_density; % mol/m^3, total mol density of polymer
g = n'*density*kB*T*NA;

Ltrue = 2*l*M/MW_monomer;
R0 = sqrt(2*l^2.*M/MW_monomer*cinf);
Rg = R0/sqrt(6);
Nk = Ltrue.^2./(2*cinf*l^2.*M/MW_monomer);
bk = Ltrue/Nk;
Nk = Nk';
L = Ltrue'./R0';
cstar = M/(NA*(2*Rg)^3)

jvals = 1:Ns;
%     hstar = 0;
sigma = -1.40*hstar^0.78;
b = 1-1.66*hstar^0.78;
intrinsic_visc = K*M.^alpha;
lambda_m = 1/S*M.*intrinsic_visc*eta_s/(NA*kB*T);
%     lambda_m = M.^1.5.*intrinsic_visc*eta_s/(NA*kB*T);
lambda_ms = lambda_m'./jvals.^(2+sigma);

% Ohnesorge number
eta_l = eta_s + sum(lambda_ms.*g);
Oh = eta_l/sqrt(rho*chi*a0)
Ec = chi/(g*a0)
beta = eta_s/eta_l

%% linear viscoelastic
nofreqs = 100;
freq = logspace(-1,3,nofreqs); % 1/s
lamlin = reshape(lambda_ms, [Nm*Ns,1]);
glin = repmat(g,[Ns,1]);
Gp = freq.^2.*glin.*lamlin.^2./(1+(lamlin.*freq).^2);
Gpp = freq.*glin.*lamlin./(1+(lamlin.*freq).^2);

Gp_PS7 = Gp;
Gpp_PS7 = Gpp;

f1 = figure(1);
pos = get(f1, 'Position');
set(f1, 'Position', [pos(1), pos(2),784.8,629.6]);
hold on
plot(freq, sum(Gp,1), '-k', 'DisplayName',"$G'_p$", 'LineWidth',3)
plot(freq, sum(Gpp,1), '--k', 'DisplayName',"$G''_p$", 'LineWidth',3)

% Gp = reshape(Gp, [Nm,Ns,nofreqs]);
% Gpp = reshape(Gpp, [Nm,Ns,nofreqs]);
% 
% for ii=1:Nm
%     plot(freq, squeeze(sum(Gp(ii,:,:),2)), ...
%         '-', 'Color', colours_MWDist(ii,:), ...
%         'HandleVisibility','off', 'LineWidth',1)
% %     plot(freq, Gp(:,ii), '-', 'Color', colours(ii,:), ...
% %         'LineWidth',1)
%     plot(freq, squeeze(sum(Gpp(ii,:,:),2)), ...
%         '--', 'Color', colours_MWDist(ii,:), ...
%         'HandleVisibility','off', 'LineWidth',1)
% end

axes1 = gca;
axes1.XScale = 'log';
axes1.YScale = 'log';

xlabel('$\omega$ (1/s)')
ylabel("$G'_p, G''_p$ (Pa)")
legend

% opts = detectImportOptions('SAOS_1000PS7.xlsx');
% SAOSData = readmatrix('SAOS_1000PS7.xlsx', opts);
opts = detectImportOptions('PS7_1000ppm_SAOS.xlsx');
SAOSData = readmatrix('PS7_1000ppm_SAOS.xlsx', opts);
freq_exp = SAOSData(:,1);
Gp_exp = SAOSData(:,2);
Gp_exp_std = SAOSData(:,3);
Gpp_exp = SAOSData(:,4) - eta_s*freq_exp;
Gpp_exp_std = SAOSData(:,5);

% complex_visc = sqrt(Gp_exp.^2 + Gpp_exp.^2)./freq_exp;

errorbar(freq_exp, Gp_exp, Gp_exp_std, 'or', 'DisplayName',"$G'$ exp", ...
    'MarkerFaceColor','r')
errorbar(freq_exp, Gpp_exp, Gpp_exp_std, 'or', 'DisplayName',"$G'' - \eta_s \omega$ exp")

omtau = freq.*lambda_m';
% Gp_ana = kB*T*density/(bk^3*Nk)...
%     *omtau.*sin(1/3*atan(omtau))./(1+omtau.^2).^(1/6);
% Gpp_ana = kB*T*density/(bk^3*Nk)...
%     *omtau.*cos(1/3*atan(omtau))./(1+omtau.^2).^(1/6);
Gp_ana = sum(g...
    .*omtau.*sin(1/3*atan(omtau))./(1+omtau.^2).^(1/6), 1);
Gpp_ana = sum(g...
    .*omtau.*cos(1/3*atan(omtau))./(1+omtau.^2).^(1/6), 1);
% plot(freq, Gp_ana, '-b', 'DisplayName',"$G'$ approx", 'LineWidth',3)
% plot(freq, Gpp_ana, '--b', 'DisplayName',"$G''$ approx", 'LineWidth',3)

% Create textbox
annotation('textbox',...
    [0.34,0.83,0.3,0.07],...
    'String',{'$c_\mathrm{PS7} = 1000$ [wppm]'},...
    'FitBoxToText','off');

xlim([5e-1, 5e2]);
ylim([1e-5, 1e1]);

%% 200 ppm PS16, no PS7

PS7_conc = 0;
PS16_conc = 200;

% project data onto same x-values
allPS_data = [PS7MW;PS16MW];
noPoints = 10000;
x = linspace(min(allPS_data), max(allPS_data),noPoints);
method = 'pchip';
yPS7 = interp1(PS7MW,PS7NormMass,x, method, 0);
yPS16 = interp1(PS16MW,PS16NormMass,x, method, 0);

mass_polymer = PS16_conc;

weight_fraction_PS16 = PS16_conc/mass_polymer;
y = yPS7*(1-weight_fraction_PS16)+yPS16*weight_fraction_PS16;
density_mass_PS = mass_polymer/vol_solvent;
maxMW = max(x);
[n, M, mol_polymer_density] = ...
    splitMWDist(graphs_on, colours_MWDist,...
    x, y, noPoints, maxMW, Nm, nu, density_mass_PS);

% get densities etc
density = mol_polymer_density; % mol/m^3, total mol density of polymer
g = n'*density*kB*T*NA;

Ltrue = 2*l*M/MW_monomer;
R0 = sqrt(2*l^2.*M/MW_monomer*cinf);
Rg = R0/sqrt(6);
Nk = Ltrue.^2./(2*cinf*l^2.*M/MW_monomer);
bk = Ltrue/Nk;
Nk = Nk';
L = Ltrue'./R0';

jvals = 1:Ns;
%     hstar = 0;
sigma = -1.40*hstar^0.78;
b = 1-1.66*hstar^0.78;
intrinsic_visc = K*M.^alpha;
lambda_m = 1/S*M.*intrinsic_visc*eta_s/(NA*kB*T);
%     lambda_m = M.^1.5.*intrinsic_visc*eta_s/(NA*kB*T);
lambda_ms = lambda_m'./jvals.^(2+sigma);

% Ohnesorge number
eta_l = eta_s + sum(lambda_ms.*g);
Oh = eta_l/sqrt(rho*chi*a0)
Ec = chi/(g*a0)
beta = eta_s/eta_l

%% linear viscoelastic
nofreqs = 100;
freq = logspace(-1,3,nofreqs); % 1/s
lamlin = reshape(lambda_ms, [Nm*Ns,1]);
glin = repmat(g,[Ns,1]);
Gp = freq.^2.*glin.*lamlin.^2./(1+(lamlin.*freq).^2);
Gpp = freq.*glin.*lamlin./(1+(lamlin.*freq).^2);

Gp_PS16 = Gp;
Gpp_PS16 = Gpp;

f2 = figure(2);
pos = get(f2, 'Position');
set(f2, 'Position', [pos(1), pos(2),784.8,629.6]);
hold on
plot(freq, sum(Gp,1), '-k', 'DisplayName',"$G'_p$", 'LineWidth',3)
plot(freq, sum(Gpp,1), '--k', 'DisplayName',"$G''_p$", 'LineWidth',3)

% Gp = reshape(Gp, [Nm,Ns,nofreqs]);
% Gpp = reshape(Gpp, [Nm,Ns,nofreqs]);
% 
% for ii=1:Nm
%     plot(freq, squeeze(sum(Gp(ii,:,:),2)), ...
%         '-', 'Color', colours_MWDist(ii,:), ...
%         'HandleVisibility','off', 'LineWidth',1)
% %     plot(freq, Gp(:,ii), '-', 'Color', colours(ii,:), ...
% %         'LineWidth',1)
%     plot(freq, squeeze(sum(Gpp(ii,:,:),2)), ...
%         '--', 'Color', colours_MWDist(ii,:), ...
%         'HandleVisibility','off', 'LineWidth',1)
% end

axes1 = gca;
axes1.XScale = 'log';
axes1.YScale = 'log';

xlabel('$\omega$ (1/s)')
ylabel("$G'_p, G''_p$ (Pa)")
legend

% opts = detectImportOptions('SAOS_200PS16.xlsx');
% SAOSData = readmatrix('SAOS_200PS16.xlsx', opts);
opts = detectImportOptions('PS16_200ppm_SAOS.xlsx');
SAOSData = readmatrix('PS16_200ppm_SAOS.xlsx', opts);
freq_exp = SAOSData(:,1);
Gp_exp = SAOSData(:,2);
Gp_exp_std = SAOSData(:,3);
Gpp_exp = SAOSData(:,4) - eta_s*freq_exp;
Gpp_exp_std = SAOSData(:,5);

% complex_visc = sqrt(Gp_exp.^2 + Gpp_exp.^2)./freq_exp;

errorbar(freq_exp, Gp_exp, Gp_exp_std, 'or', 'DisplayName',"$G'$ exp", ...
    'MarkerFaceColor','r')
errorbar(freq_exp, Gpp_exp, Gpp_exp_std, 'or', 'DisplayName',"$G'' - \eta_s \omega$ exp")


% omtau = freq*lambda_m(1);
omtau = freq.*lambda_m';
% Gp_ana = kB*T*density/(bk^3*Nk)...
%     *omtau.*sin(1/3*atan(omtau))./(1+omtau.^2).^(1/6);
% Gpp_ana = kB*T*density/(bk^3*Nk)...
%     *omtau.*cos(1/3*atan(omtau))./(1+omtau.^2).^(1/6);
Gp_ana = sum(g...
    .*omtau.*sin(1/3*atan(omtau))./(1+omtau.^2).^(1/6), 1);
Gpp_ana = sum(g...
    .*omtau.*cos(1/3*atan(omtau))./(1+omtau.^2).^(1/6), 1);
% plot(freq, Gp_ana, '-b', 'DisplayName',"$G'$ approx", 'LineWidth',3)
% plot(freq, Gpp_ana, '--b', 'DisplayName',"$G''$ approx", 'LineWidth',3)

% Create textbox
annotation('textbox',...
    [0.34,0.83,0.3,0.07],...
    'String',{'$c_\mathrm{PS16} = 200$ [wppm]'},...
    'FitBoxToText','off');

xlim([5e-1, 5e2]);
ylim([1e-5, 1e1]);

%% plot both together
%%%%%%%%%%%%%%%%%%%%%% COMBINE AND PLOT
Gp = Gp_PS7 + Gp_PS16;
Gpp = Gpp_PS7 + Gpp_PS16;

f3 = figure(3);
pos = get(f3, 'Position');
set(f3, 'Position', [pos(1), pos(2),784.8,629.6]);
hold on
plot(freq, sum(Gp,1), '-k', 'DisplayName',"$G'$", 'LineWidth',3)
plot(freq, sum(Gpp,1), '--k', 'DisplayName',"$G''$", 'LineWidth',3)

% Gp = reshape(Gp, [Nm,Ns,nofreqs]);
% Gpp = reshape(Gpp, [Nm,Ns,nofreqs]);
% 
% for ii=1:Nm
%     plot(freq, squeeze(sum(Gp(ii,:,:),2)), ...
%         '-', 'Color', colours_MWDist(ii,:), ...
%         'HandleVisibility','off', 'LineWidth',1)
% %     plot(freq, Gp(:,ii), '-', 'Color', colours(ii,:), ...
% %         'LineWidth',1)
%     plot(freq, squeeze(sum(Gpp(ii,:,:),2)), ...
%         '--', 'Color', colours_MWDist(ii,:), ...
%         'HandleVisibility','off', 'LineWidth',1)
% end

axes1 = gca;
axes1.XScale = 'log';
axes1.YScale = 'log';

xlabel('$\omega$ (1/s)')
ylabel("$G', G''$ (Pa)")
legend

opts = detectImportOptions('SAOS_mix_1000PS7_200PS16.csv');
SAOSData = readmatrix('SAOS_mix_1000PS7_200PS16.csv', opts);
freq_exp = SAOSData(:,1);
Gp_exp = SAOSData(:,2);
Gpp_exp = SAOSData(:,3) - eta_s*freq_exp;

% complex_visc = sqrt(Gp_exp.^2 + Gpp_exp.^2)./freq_exp;

% errorbar(freq_exp, Gp_exp, Gp_exp_std, 'or', 'DisplayName',"$G'$ exp", ...
%     'MarkerFaceColor','r')
% errorbar(freq_exp, Gpp_exp, Gpp_exp_std, 'or', 'DisplayName',"$G''$ exp")

plot(freq_exp, Gp_exp, 'or', 'DisplayName',"$G'$ exp", ...
    'MarkerFaceColor','r')
plot(freq_exp, Gpp_exp, 'or', 'DisplayName',"$G''$ exp")


% omtau = freq*lambda_m(1);
omtau = freq.*lambda_m';
% Gp_ana = kB*T*density/(bk^3*Nk)...
%     *omtau.*sin(1/3*atan(omtau))./(1+omtau.^2).^(1/6);
% Gpp_ana = kB*T*density/(bk^3*Nk)...
%     *omtau.*cos(1/3*atan(omtau))./(1+omtau.^2).^(1/6);
Gp_ana = sum(g...
    .*omtau.*sin(1/3*atan(omtau))./(1+omtau.^2).^(1/6), 1);
Gpp_ana = sum(g...
    .*omtau.*cos(1/3*atan(omtau))./(1+omtau.^2).^(1/6), 1);
% plot(freq, Gp_ana, '-b', 'DisplayName',"$G'$ approx", 'LineWidth',3)
% plot(freq, Gpp_ana, '--b', 'DisplayName',"$G''$ approx", 'LineWidth',3)

% Create textbox
annotation('textbox',...
    [0.34,0.83,0.3,0.07],...
    'String',{'$c_\mathrm{PS16} = 200$ [wppm]', ...
              '$c_\mathrm{PS7} = 1000$ [wppm]'},...
    'FitBoxToText','off');

% xlim([5e-1, 5e2]);
% ylim([1e-5, 1e1]);

% sometimes figures are out of frame
findfigs














