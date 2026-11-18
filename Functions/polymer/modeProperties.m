function pm = modeProperties(modes, p, Ns, varargin)
%MODEPROPERTIES Moduli, extensibilities and relaxation times of each mode.
%
%   pm = modeProperties(modes, p, Ns, Name, Value) turns the modes from
%   twoSpeciesModes/blendModes (fields n, M, molDensity) into the
%   multi-mode FENE-PM input of solveFilamentThinning:
%       g_i       = n_i kB T                                     (Eq. 15)
%       L_i^2     = 3 (R_max^2/<R^2>_eq - 1)                     (Eqs. 23, 24)
%       N_K,i     = (L_i^2 + 1)/3   Kuhn steps, for the CDD model (Eq. 28)
%       tau_Z,i   = M_i [eta] eta_s / (U_tau_eta N_A kB T)       (Eq. 17)
%       tau_is    = tau_Z,i / s^(2+p),  p = -1.40 h*^0.78        (Eqs. 19, 20)
%
%   Options (defaults in brackets)
%     'hstar'  [p.hstar]  HI parameter of the Zimm spectrum
%     'S'      [p.S]      U_tau_eta; S = 1 gives the 2.36x longer times the
%                         original Figs 3, 8 and 9 used
%
%   Fields (column vectors are Nm-by-1): g [Pa], L, Nk, lambda_ms (Nm-by-Ns),
%   lambda_m (tau_Z, 1-by-Nm), Ltrue (R_max), R0 (sqrt(<R^2>_eq)) (1-by-Nm).

opt = parseOptions(struct('hstar', p.hstar, 'S', p.S), varargin{:});

M = modes.M;
pm.g = modes.n'*modes.molDensity*p.kB*p.T*p.NA;

pm.Ltrue = 2*p.l*M/p.MW_monomer;
pm.R0 = sqrt(2*p.l^2.*M/p.MW_monomer*p.cinf);
pm.L = sqrt(3*((pm.Ltrue'./pm.R0').^2-1));
pm.Nk = (pm.L.^2 + 3)/3;

jvals = 1:Ns;
sigma = -1.40*opt.hstar^0.78;
intrinsic_visc = p.K*M.^p.alpha;
pm.lambda_m = 1/opt.S*M.*intrinsic_visc*p.eta_s/(p.NA*p.kB*p.T);
pm.lambda_ms = pm.lambda_m'./jvals.^(2+sigma);
end
