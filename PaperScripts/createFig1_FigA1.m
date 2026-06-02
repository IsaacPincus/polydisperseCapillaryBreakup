clear variables;
close all;
warning('off','MATLAB:nearlySingularMatrix');

datapath = "../Data";
addpath(datapath)
addpath("../Functions")

% PS constants
alpha = 0.5;
K = 8e-8; % m^3/g, Mark-Houwink paramceter, fitted
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
rho = 1000; % kg/m^3, fluid density
alpha_C1 = 0.4;
C1 = 9/(4*alpha_C1^3);
rho = rho*C1*8;
tmax = 100; % s
cutoff = 1e-4;
Ns = 20;
PS7_conc = 1000;
PS16_conc = 500;
% PS16_conc = 200;
hstar = 0.2671;
mu = eta_s;
a0 = 4.8e-4;
% a0 = 20e-4;
X = 1;
colours_MWDist = parula(6);
colours_MWDist = colours_MWDist([1,5],:);
graphs_on = true;

filePS7 = fullfile(datapath, "PS7MWDist.csv");
opts = detectImportOptions(filePS7);
dataPS7 = readmatrix(filePS7, opts);
PS7MW = dataPS7(:,1);
PS7NormMass = dataPS7(:,2);

filePS16 = fullfile(datapath, "PS16MWDist.csv");
opts = detectImportOptions(filePS16);
dataPS16 = readmatrix(filePS16, opts);
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

figure(2);
hold on
copyobj(figure(2).Children.Children, figure(1).Children);

close 2

f1 = gcf;
objs = f1.Children.Children;

objs(1).Color = colours_MWDist(1,:);
objs(2).FaceColor = colours_MWDist(1,:);
% objs(3).Color = colours_MWDist(1,:);
% objs(3).LineWidth = 3;

objs(4).Color = colours_MWDist(2,:);
objs(5).FaceColor = colours_MWDist(2,:);
% objs(6).Color = colours_MWDist(2,:);
% objs(6).LineWidth = 3;

Nm = 2;
mol_polymer_density = mol_polymer_densityPS7 + mol_polymer_densityPS16;
mol_frac_PS7 = mol_polymer_densityPS7/mol_polymer_density;
n = [nPS16*(1-mol_frac_PS7), nPS7*mol_frac_PS7];
M = [MPS16, MPS7];

objs(1).YData = [0,objs(4).YData(2)*n(1)/n(2)];

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

% Create textarrow
annotation(gcf,'textbox',[0.361660079051385,0.782152230971129, ...
    0.0800395256917,0.05748031496063],'String',{'PS7'});
annotation(gcf,'textarrow',[0.7549,0.59930],...
    [0.28740,0.21378],'String',{'$M_{L,i=1}$'})

% Create textarrow
annotation(gcf,'textbox',[0.65316205533597,0.568241469816274, ...
    0.0899209486166,0.058792650918636],'String',{'PS16'})
annotation(gcf,'textarrow',[0.4457,0.3454],...
    [0.60412,0.54055],'String',{'$M_{L,i=2}$'})

% Create textbox
annotation(gcf,'textbox',...
    [0.586956,0.769028,0.27470355,0.11653543],...
    'String',{'$c_\mathrm{PS7} = 1000$ [wppm]', ...
    '$c_\mathrm{PS16} = 500$ [wppm]'});

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables;

% PS constants
alpha = 0.5;
K = 11e-8; % m^3/g, Mark-Houwink parameter, fitted
% alpha = 1;
% K = 5e-11; % m^3/g, Mark-Houwink parameter, fitted
density_solvent = 0.98*997e3; % g/m^3
vol_solvent = 1e6/density_solvent; % vol in m^3 of 1Mil g of solvent
eta_s = 0.056; % Pa.s, solvent viscosity
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
alpha_C1 = 0.3;
C1 = 9/(4*alpha_C1^3);
rho = rho*C1*8;
tmax = 100; % s
cutoff = 1e-4;
Nm = 30;
Ns = 1;
PS7_conc = 1000;
PS16_conc = 500;
hstar = 0.2671;
mu = eta_s;
a0 = 4.8e-4;
X = 1;
colours_MWDist = parula(Nm);
graphs_on = true;

opts = detectImportOptions("PS7MWDist.csv");
dataPS7 = readmatrix("PS7MWDist.csv", opts);
PS7MW = dataPS7(:,1);
PS7NormMass = dataPS7(:,2);

opts = detectImportOptions("PS16MWDist.csv");
dataPS16 = readmatrix("PS16MWDist.csv", opts);
PS16MW = dataPS16(:,1);
PS16NormMass = dataPS16(:,2);

% %% fit together, multiple modes

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

% Create textbox
annotation(gcf,'textbox',...
    [0.5 0.75 0.4 0.15],...
    'String',{'$c_\mathrm{PS7} = 1000$ [wppm]', ...
    '$c_\mathrm{PS16} = 500$ [wppm]','$N_i = 30$'});

% Create textarrow
annotation(gcf,'textarrow',[0.472332015810277 0.350790513833992],...
    [0.583224115334207 0.516382699868938],'String',{'number density'});

% Create textarrow
annotation(gcf,'textarrow',[0.563241106719368 0.596837944664032],...
    [0.397116644823067 0.292267365661861],'String',{'mass density'})







































