function p = defaultParameters()
%DEFAULTPARAMETERS Shared physical constants and PS/DOP model parameters.
%
%   p = defaultParameters() returns the parameters of Table I of the paper,
%   used by every figure. Only numerical settings (the number of Zimm
%   modes Ns, concentrations, sweep ranges) are chosen per figure, at the
%   top of the corresponding computeFig* function.

% --- Physical constants ---------------------------------------------------
p.NA = 6.02e23;          % 1/mol, Avogadro's constant
p.kB = 1.38e-23;         % J/K, Boltzmann constant
p.T  = 298;              % K, temperature

% --- Polystyrene ----------------------------------------------------------
p.l = 1.54e-10;          % m, C-C bond length
p.cinf = 9.7;            % characteristic ratio
p.MW_monomer = 104;      % g/mol, monomer molecular weight
p.K = 8e-8;              % m^3/g, Mark-Houwink prefactor (fitted)
p.alpha = 0.5;           % Mark-Houwink exponent
p.nu = 0.5;              % exponent of the extensibility-averaged MW (splitMWDist)

% --- Solvent (DOP) and filament -------------------------------------------
p.eta_s = 0.059;                       % Pa.s, solvent viscosity
p.density_solvent = 0.98*997e3;        % g/m^3
p.vol_solvent = 1e6/p.density_solvent; % m^3 of 1e6 g of solvent (conc. in wppm)
p.chi = 3.1e-2;          % N/m, surface tension
p.X = 1;                 % capillary pressure is (2X-1)*chi/a
p.rhoFluid = 985;        % kg/m^3, fluid density
p.rho = inertialDensity(p.rhoFluid, 0.25); % density scaled by the inertia constant C = 144
p.a0 = 6.2e-4;           % m, initial filament radius
p.tmax = 100;            % s, maximum simulated time
p.cutoff = 1e-3;         % simulation stops when a/a0 reaches this value

% --- Polymer relaxation spectrum and drag ---------------------------------
p.hstar = 0.26;          % hydrodynamic interaction parameter of the Zimm spectrum
p.hsk = 0.027;           % Kuhn-segment HI parameter used by zetacal (CDD)
p.phieq = 0;             % c/c* used by zetacal (0: no concentration dependence)
p.S = 2.36;              % U_tau_eta in the Zimm time (Eq. 17 of the paper)

% --- Molecular weight distribution grid -----------------------------------
p.noPoints = 10000;      % points on the common MW grid
p.interpMethod = 'pchip';
end
