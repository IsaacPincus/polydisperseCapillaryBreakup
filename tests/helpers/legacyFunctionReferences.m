function [cases, unusedCases] = legacyFunctionReferences(legacyDir)
%LEGACYFUNCTIONREFERENCES Inputs and ORIGINAL outputs for the function tests.
%
%   Must run with only the original code on the path (generateReferenceData
%   arranges this). Each case is a struct with fields inputs and outputs.
%   Where a refactored function replaced code that sat inside the original
%   paper scripts, that code is reproduced below verbatim and marked
%   "from createFig*.m". Cases are chosen where the paper-consistent code
%   should still agree with the original (see the notes at each case).
%
%   unusedCases covers the functions no paper script calls.

cases = struct();
legacyData = fullfile(legacyDir, 'Data');

% --- shared constants, as in the original scripts -----------------------------
NA = 6.02e23; kB = 1.38e-23; T = 298; l = 1.54e-10; cinf = 9.7;
MW_monomer = 104; K = 8e-8; alpha = 0.5; nu = 0.5; eta_s = 0.059;
density_solvent = 0.98*997e3; vol_solvent = 1e6/density_solvent;
chi = 3.1e-2; X = 1; mu = eta_s;
rho = 1000; alpha_C1 = 0.25; C1 = 9/(4*alpha_C1^3); rho = rho*C1*8;
tmax = 100; cutoff = 1e-3; a0 = 6.3e-4; hsk = 0.027; phieq = 0;

% --- inertialDensity (from createFig4_Fig5.m) ---------------------------------
rho985 = 985; C1 = 9/(4*0.25^3); rho985 = rho985*C1*8;
cases.inertialDensity = io(struct('rho', 985, 'alpha_C1', 0.25), struct('rhoEff', rho985));

% --- ubcal ----------------------------------------------------------------------
in = struct('hstar', {{0.027, 0.2671, 0}}, 'nkperb', [1, 2.5, 10, 1e3, 3.3e4]);
out.Ub = cellfun(@(h) ubcal(h, in.nkperb), in.hstar, 'UniformOutput', false);
cases.ubcal = io(in, out);

% --- zetacal: every model on a grid of stretch and concentration ---------------
[Qg, Pg, Hg, Mg] = ndgrid([1, 2.5, 20, 150, 400], [0, 0.3, 3, 1e3], [0.027, 0], 0:3);
in = struct('Q', Qg(:), 'phieq', Pg(:), 'hsk', Hg(:), 'himodel', Mg(:), 'nk', 1e5);
out = struct('zet', nan(size(in.Q)), 'G', nan(size(in.Q)), 'zb', nan(size(in.Q)), 'Nb', nan(size(in.Q)));
for k = 1:numel(in.Q)
    [out.zet(k), G, zb, Nb] = zetacal(in.Q(k), 1, in.phieq(k), in.nk, in.hsk(k), in.himodel(k));
    out.G(k) = G; out.zb(k) = zb; out.Nb(k) = Nb;
end
cases.zetacal = io(in, out);

% --- loadMWDistributions (from createFig3.m) ------------------------------------
opts = detectImportOptions(fullfile(legacyData, 'PS7MWDist.csv'));
dataPS7 = readmatrix(fullfile(legacyData, 'PS7MWDist.csv'), opts);
PS7MW = dataPS7(:,1); PS7NormMass = dataPS7(:,2);
opts = detectImportOptions(fullfile(legacyData, 'PS16MWDist.csv'));
dataPS16 = readmatrix(fullfile(legacyData, 'PS16MWDist.csv'), opts);
PS16MW = dataPS16(:,1); PS16NormMass = dataPS16(:,2);
allPS_data = [PS7MW;PS16MW];
noPoints = 10000;
x = linspace(min(allPS_data), max(allPS_data),noPoints);
method = 'pchip';
yPS7 = interp1(PS7MW,PS7NormMass,x, method, 0);
yPS16 = interp1(PS16MW,PS16NormMass,x, method, 0);
cases.loadMWDistributions = io(struct(), struct('x', x, 'yPS7', yPS7, 'yPS16', yPS16));

% --- splitMWDist ------------------------------------------------------------------
colours = parula(30);
in = struct('N', {{1, 5, 30}}, 'wPS16', 0.3, 'density', 1300/vol_solvent);
in.x = x;
in.y = yPS7*(1-in.wPS16) + yPS16*in.wPS16;
for k = 1:numel(in.N)
    [n, M, mol] = splitMWDist(false, colours, in.x, in.y, noPoints, max(x), in.N{k}, nu, in.density);
    out = struct('n', n, 'M', M, 'molDensity', mol);
    outs{k} = out; %#ok<AGROW>
end
cases.splitMWDist = io(in, struct('results', {outs}));

% --- twoSpeciesModes (from createFig3.m / get_outputs) -------------------------
conc = [1000 500; 200 0; 0 200; 1000 100];
outs = cell(1, size(conc, 1));
for k = 1:size(conc, 1)
    outs{k} = legacyTwoSpecies(x, yPS7, yPS16, conc(k,1), conc(k,2), noPoints, nu, vol_solvent);
end
cases.twoSpeciesModes = io(struct('conc', conc), struct('results', {outs}));

% --- blendModes (from createFig7.m) ---------------------------------------------
conc = [1000 100 30; 1000 500 30; 200 300 80; 1000 0 1];
outs = cell(1, size(conc, 1));
for k = 1:size(conc, 1)
    PS7_conc = conc(k,1); PS16_conc = conc(k,2); Nm = conc(k,3);
    mass_polymer = PS16_conc + PS7_conc;
    weight_fraction_PS16 = PS16_conc/mass_polymer;
    y = yPS7*(1-weight_fraction_PS16)+yPS16*weight_fraction_PS16;
    density_mass_PS = mass_polymer/vol_solvent;
    [n, M, mol_polymer_density] = splitMWDist(false, parula(Nm), x, y, noPoints, max(x), Nm, nu, density_mass_PS);
    outs{k} = struct('n', n, 'M', M, 'molDensity', mol_polymer_density, 'y', y);
end
cases.blendModes = io(struct('conc', conc), struct('results', {outs}));

% --- modeProperties (from createFig2.m, createFig3.m, get_outputs, createFig8.m)
modes = legacyTwoSpecies(x, yPS7, yPS16, 1000, 100, noPoints, nu, vol_solvent);
variants = {  % Ns, hstar, S
    30, 0.267, 1
    10, 0.2671, 2.36
    50, 0.26, 2.36
    3, 0.26, 2.36};
outs = cell(1, size(variants, 1));
for k = 1:size(variants, 1)
    [Ns, hstar, S] = variants{k,:};
    outs{k} = legacyModeProperties(modes.n, modes.M, modes.molDensity, Ns, hstar, S, 1, 'FENE-P', ...
        NA, kB, T, l, cinf, MW_monomer, K, alpha, eta_s);
    % Nk is now (L^2+1)/3 as in the paper, not Lmax^2/R0^2
    outs{k} = rmfield(outs{k}, 'Nk');
end
cases.modeProperties = io(struct('modes', modes, 'variants', {variants}), struct('results', {outs}));

% --- linearViscoelasticModuli (from createFig2.m) --------------------------------
pm = outs{3};
freq = logspace(-1,3,100);
Nm = numel(pm.g); Ns = size(pm.lambda_ms, 2);
lamlin = reshape(pm.lambda_ms, [Nm*Ns,1]);
glin = repmat(pm.g,[Ns,1]);
Gp = freq.^2.*glin.*lamlin.^2./(1+(lamlin.*freq).^2);
Gpp = freq.*glin.*lamlin./(1+(lamlin.*freq).^2);
cases.linearViscoelasticModuli = io(struct('pm', pm, 'freq', freq), struct('Gp', Gp, 'Gpp', Gpp));

% --- pack/unpackConformation (loops from the original solvers) ----------------
Nm = 3; Ns = 4;
v = (1:2*Nm*Ns)'/7;
count = 0;
for mm=1:Nm
    for ss = 1:Ns
        count = count + 1;
        Az(mm,ss) = v(count);
        count = count + 1;
        Ar(mm,ss) = v(count);
    end
end
cases.packConformation = io(struct('v', v, 'Nm', Nm, 'Ns', Ns), struct('Az', Az, 'Ar', Ar));

% --- solver cases -------------------------------------------------------------------
% 1000 wppm PS7 + 100 wppm PS16. The refactored solver follows the paper's
% Eq. 22 (f = L^2/(L^2 + 3Ns - sum tr A)) and uses the FENE trace Az + 2Ar
% throughout; the original solvers used L^2 + 3 and, with conformation-
% dependent drag, Az - 2Ar. They agree for constant drag with Ns = 1, which
% is what is compared here (SolverPhysicsTest covers the rest).
small = legacyModeProperties(modes.n, modes.M, modes.molDensity, 3, 0.2671, 2.36, 1, 'FENE-P', ...
    NA, kB, T, l, cinf, MW_monomer, K, alpha, eta_s);
oneMode = legacyModeProperties(modes.n, modes.M, modes.molDensity, 1, 0.26, 2.36, 1, 'FENE-P', ...
    NA, kB, T, l, cinf, MW_monomer, K, alpha, eta_s);
p = struct('chi', chi, 'X', X, 'eta_s', mu, 'rho', rho, 'a0', a0, ...
    'cutoff', cutoff, 'tmax', tmax, 'hsk', hsk, 'phieq', phieq);
names = {'a','Az','Ar','f','stress','ep','Wi','t','te','ye','ie'};

single = oneMode;  % PS7 only
single.g = oneMode.g(2); single.L = oneMode.L(2); single.lambda_ms = oneMode.lambda_ms(2,:); single.Nk = oneMode.Nk(2);

out = cell(1, numel(names));
[out{:}] = filamentThinningFENEPM_inertia(chi, X, mu, rho, 2, 1, oneMode.L, oneMode.lambda_ms, oneMode.g, a0, cutoff, tmax);
cases.solverConstantDrag = io(struct('p', p, 'modes', oneMode), cell2struct(out, names, 2));

[out{:}] = filamentThinningFENEPM_inertia(chi, X, mu, rho, 1, 1, single.L, single.lambda_ms, single.g, a0, cutoff, tmax);
cases.solverConstantDragSingleMode = io(struct('p', p, 'modes', single), cell2struct(out, names, 2));

a0i = 4e-4;
y0 = [a0i; repmat([3; 0.6], 2, 1)];
epInitial = 30;
pInitial = p; pInitial.a0 = a0i;
[out{:}] = filamentThinningFENEPM_inertia_initial(chi, X, mu, rho, 2, 1, oneMode.L, oneMode.lambda_ms, oneMode.g, a0i, cutoff, tmax, y0, epInitial);
cases.solverInitialState = io(struct('p', pInitial, 'modes', oneMode, 'conformation', y0(2:end), ...
    'epInitial', epInitial), cell2struct(out, names, 2));

% --- simulatedRelaxationTime (from createFig3/4/8.m and Fig9_*.m) --------------
s = cases.solverConstantDrag.outputs;
t = s.t; a = s.a; ep = s.ep;
outs = struct();
% get_outputs (Figs 4-6): threshold 1.2, unshifted fit
t_max_ep_all = t(islocalmax(ep)==1);
t_c = t_max_ep_all(1);
ep_min_array = islocalmin(ep)==1;
ep_min_all = ep(ep_min_array&(t>t_c));
ep_min = min(ep_min_all);
tau_ep_min = 2/(3*ep_min);
indices_EC = find((ep<1.2*ep_min)&(t>t_c));
[xxx,yyy] = prepareCurveData(t(indices_EC),a(indices_EC));
fitobject = fit(xxx, yyy, 'exp1');
outs.getOutputs = struct('t_c', t_c, 'ep_min', ep_min, 'tau_ep_min', tau_ep_min, ...
    'indices_EC', indices_EC, 'tau_EC', -1./(3*fitobject.b));
% createFig3: shifted fit, t_ep_min
t_ep_min = min(t(ep_min_array&(t>t_c)));
[xxx,yyy] = prepareCurveData(t(indices_EC)-t_c,a(indices_EC));
fitobject = fit(xxx, yyy, 'exp1');
outs.shifted = struct('t_c', t_c, 't_ep_min', t_ep_min, 'tau_EC', -1./(3*fitobject.b), ...
    'fitA', fitobject.a);
% createFig8: threshold 2, index of the minimum
indices_ep_min = find(ep_min_array&(t>t_c));
[ep_min, index] = min(ep(indices_ep_min));
ep_min_index = indices_ep_min(index);
indices_EC = find((ep<2*ep_min)&(t>t_c));
[xxx,yyy] = prepareCurveData(t(indices_EC),a(indices_EC));
fitobject = fit(xxx, yyy, 'exp1');
outs.threshold2 = struct('ep_min_index', ep_min_index, 'indices_EC', indices_EC, ...
    'tau_EC', -1./(3*fitobject.b));
% Fig9_FastRetraction_CDD: t_c = 0, global minimum
t_c = 0;
[ep_min, index_ep_min] = min(ep);
t_ep_min = t(index_ep_min);
indices_EC = find((ep<1.2*ep_min)&(t>t_c));
[xxx,yyy] = prepareCurveData(t(indices_EC)-t_c,a(indices_EC));
fitobject = fit(xxx, yyy, 'exp1');
outs.globalMinimum = struct('ep_min', ep_min, 't_ep_min', t_ep_min, ...
    'tau_ep_min', 2/(3*ep_min), 'tau_EC', -1./(3*fitobject.b));
cases.simulatedRelaxationTime = io(struct('t', t, 'a', a, 'ep', ep), outs);

% --- readThinningCurves and experimentalRelaxationTime (from createFig4_Fig5.m,
%     createFig6a.m and createFig6b.m) --------------------------------------------
spreadsheetName = fullfile(legacyData, 'ThinningCurves_PS7_PS16_andMixtures.xlsx');
opts = detectImportOptions(spreadsheetName, 'Sheet','PS16');
opts.SelectedVariableNames = [1:2,5:6,9:10,13:14,17:18,21:22];
M = readmatrix(spreadsheetName,opts);
t16 = M(2:end,1:2:11);
D16 = M(2:end,2:2:12);
opts = detectImportOptions(spreadsheetName, 'Sheet','200ppmPS7+PS16');
opts.SelectedVariableNames = [1:2];
M = readmatrix(spreadsheetName,opts);
t7 = M(:,1); D7 = M(:,2);
opts = detectImportOptions(spreadsheetName, 'Sheet','200ppmPS7+PS16');
opts.SelectedVariableNames = [1:2,5:6,9:10,13:14,17:18,21:22,25:26];
M = readmatrix(spreadsheetName,opts);
tMix = M(2:end,1:2:13);
DMix = M(2:end,2:2:14);
cases.readThinningCurves = io(struct(), struct('t16', t16, 'D16', D16, 't7', t7, 'D7', D7, ...
    'tMix', tMix, 'DMix', DMix));

outs = struct();
% Fig 4a: radius, threshold 5, start point
te = t16(:,1); ae = D16(:,1)/2;
tavg = (te(2:end) + te(1:end-1))/2;
aavg = (ae(2:end) + ae(1:end-1))/2;
ep_exp = -2./aavg.*diff(ae)./diff(te);
[ep_min_exp, index_ep_min] = min(ep_exp(tavg>0)); %#ok<ASGLU>
tau_ep_min_exp = 2/(3*ep_min_exp);
indices_EC = find((ep_exp<5*ep_min_exp)&(tavg>0));
[xxx,yyy] = prepareCurveData(tavg(indices_EC),aavg(indices_EC));
fitobject = fit(xxx, yyy, 'exp1', 'StartPoint', [0.1,-1/(3*tau_ep_min_exp)]);
outs.radiusStartPoint = struct('tavg', tavg, 'ep', ep_exp, 'tau_ep_min', tau_ep_min_exp, ...
    'tau_EC', -1./(3*fitobject.b));
% Fig 4b: radius, threshold 3
te = t7; ae = D7/2;
tavg = (te(2:end) + te(1:end-1))/2;
aavg = (ae(2:end) + ae(1:end-1))/2;
ep_exp = -2./aavg.*diff(ae)./diff(te);
[ep_min_exp, index_ep_min] = min(ep_exp(tavg>0)); %#ok<ASGLU>
tau_ep_min_exp = 2/(3*ep_min_exp);
indices_EC = find((ep_exp<3*ep_min_exp)&(tavg>0));
[xxx,yyy] = prepareCurveData(tavg(indices_EC),aavg(indices_EC));
fitobject = fit(xxx, yyy, 'exp1');
outs.radiusThreshold3 = struct('tau_ep_min', tau_ep_min_exp, 'tau_EC', -1./(3*fitobject.b));
% Fig 6b (200 wppm PS7): diameter, positive rates only
te = tMix(:,3); De = DMix(:,3);
tavg = (te(2:end) + te(1:end-1))/2;
Davg = (De(2:end) + De(1:end-1))/2;
ep_exp = -2./Davg.*diff(De)./diff(te);
[ep_min_exp, index_ep_min] = min(ep_exp(tavg>0&ep_exp>0)); %#ok<ASGLU>
tau_ep_min_exp = 2/(3*ep_min_exp);
indices_EC = find((ep_exp<5*ep_min_exp)&(tavg>0));
[xxx,yyy] = prepareCurveData(tavg(indices_EC),Davg(indices_EC));
fitobject = fit(xxx, yyy, 'exp1');
outs.diameterPositive = struct('ep', ep_exp, 'tau_ep_min', tau_ep_min_exp, ...
    'tau_EC', -1./(3*fitobject.b));
cases.experimentalRelaxationTime = io(struct('t16', t16(:,1), 'D16', D16(:,1), 't7', t7, ...
    'D7', D7, 'tMix', tMix(:,3), 'DMix', DMix(:,3)), outs);

% --- plateRetraction (from Fig9_FastRetraction_CDD.m) ----------------------------
zeta0 = 0.5; xi0 = 1; R0 = 3e-3; L0 = 1e-3; Lend = 1.3e-3; V0 = 0.5e-3;
Hend = Lend/L0; Lam0 = L0/R0;
Hfunc = @(t) 1 + t;
Hderivfunc = @(t) 1;
tpull = fzero(@(t) Hfunc(t) - Hend, 1)*L0/V0;
options = odeset('RelTol', 1e-6,'AbsTol',1e-9);
Nm = 2; Ns = 1;   % Eq. 22 and the original agree for Ns = 1
y0 = ones(Nm*Ns*2,1);
sol = ode15s(@(t,y) pullODE(t, y, Nm, Ns, oneMode.L, oneMode.lambda_ms,...
    Hfunc, Hderivfunc, zeta0, Lam0, V0, L0), [0, tpull], y0, options);
points = 1000;
[y,yp] = deval(sol,linspace(0,tpull,points)); %#ok<ASGLU>
H = Hfunc(tpull*V0/L0);
dH = Hderivfunc(tpull*V0/L0);
z0 = ((zeta0-1)*zeta0)/(1-2*zeta0)^2;
Rpull = (xi0*((H/z0 - 4)./(H/z0 - H*4)).^(3/(4*Lam0^2)))*R0;
epend = (-3/(2*Lam0^2)*1./(H/(4*z0) - 1).*dH./H)*V0/L0;
geom = struct('R0', R0, 'L0', L0, 'Lend', Lend, 'zeta0', zeta0, 'xi0', xi0);
cases.plateRetraction = io(struct('pm', oneMode, 'V0', V0, 'geom', geom), ...
    struct('tpull', tpull, 'Rpull', Rpull, 'epend', epend, 'conformation', y(:,end)));

% =============================================================================
% Functions that no paper script calls (Functions/unused)
% =============================================================================
unusedCases = struct();
texp = t16(:,1); Dexp = D16(:,1)*1e-3;   % PS16 5 wppm, SI units

unusedCases.getInitialStretch = tryCase(struct('L0', 1e-3, 'L1', 1.3e-3, 'R0', 3e-3), ...
    @(in) namedOutputs({'Az','Ar'}, @getInitialStretch, in.L0, in.L1, in.R0));

s = cases.solverConstantDrag.outputs;
unusedCases.getTimes = tryCase(struct('sol', s, 'chi', chi, 'mu', mu, 'N', 2, 'L', oneMode.L, ...
    'lambda', oneMode.lambda_ms, 'g', oneMode.g, 'a0', a0), ...
    @(in) namedOutputs({'tau_Simon','longest_relax','longest_relax_EC','t_c','tau_EC', ...
        'longest_relax_min_ep','tau_EC_fit','tau_EC_gt_v','t_breakup'}, @getTimes, ...
        in.sol.a, in.sol.Az, in.sol.Ar, in.sol.f, in.sol.stress, in.sol.ep, in.sol.Wi, ...
        in.sol.t, in.sol.te, in.sol.ye, in.sol.ie, in.chi, in.mu, in.N, in.L, in.lambda, ...
        in.g, in.a0, false));

unusedCases.fitGaussiansToMWDist = tryCase(struct('x', x, 'y', yPS7*0.7 + yPS16*0.3, ...
    'density', 1300/vol_solvent), ...
    @(in) namedOutputs({'n','M','molDensity'}, @fitGaussiansToMWDist, false, parula(2), ...
        in.x, in.y, numel(in.x), max(in.x), 2, nu, in.density));

unusedCases.getStretchExperiments = tryCase(struct('texp', texp, 'Dexp', Dexp, ...
    'L', small.L, 'lambda', small.lambda_ms), ...
    @(in) namedOutputs({'y0','a0'}, @getStretchExperiments, in.texp, in.Dexp, 2, 3, in.L, in.lambda));

unusedCases.getStretchExperimentsStressBalanced = tryCase(struct('texp', texp, 'Dexp', Dexp, ...
    'p', p, 'modes', small), ...
    @(in) namedOutputs({'t0','y0','a0'}, @getStretchExperimentsStressBalanced, in.p.chi, ...
        in.p.X, in.p.eta_s, in.p.rho, 2, 3, in.modes.L, in.modes.lambda_ms, in.modes.g, ...
        in.texp, in.Dexp));

unusedCases.filamentThinningFENEPM_inertia_initialEXP = tryCase(struct('texp', texp, ...
    'Dexp', Dexp, 'p', p, 'modes', small), ...
    @(in) namedOutputs(names, @filamentThinningFENEPM_inertia_initialEXP, in.p.chi, in.p.X, ...
        in.p.eta_s, in.p.rho, 2, 3, in.modes.L, in.modes.lambda_ms, in.modes.g, ...
        in.p.cutoff, in.p.tmax, in.texp, in.Dexp));

unusedCases.filamentThinningFENEPM_inertia_simple_CDD = tryCase(struct('p', p, 'modes', small, ...
    'zet_rat', 2), ...
    @(in) namedOutputs(names, @filamentThinningFENEPM_inertia_simple_CDD, in.p.chi, in.p.X, ...
        in.p.eta_s, in.p.rho, 2, 3, in.modes.L, in.modes.Nk, in.zet_rat, in.modes.lambda_ms, ...
        in.modes.g, in.p.a0, in.p.cutoff, in.p.tmax));
end

% =============================================================================
function c = io(inputs, outputs)
c = struct('inputs', inputs, 'outputs', outputs);
end

function c = tryCase(inputs, run)
% unused functions may fail on any input; record the error instead
c = struct('inputs', inputs, 'outputs', [], 'error', '');
try
    c.outputs = run(inputs);
catch err
    c.error = err.message;
    fprintf('  (original function failed: %s)\n', err.message);
end
end

function s = namedOutputs(names, fun, varargin)
out = cell(1, numel(names));
[out{:}] = fun(varargin{:});
s = cell2struct(out, names, 2);
end

function m = legacyTwoSpecies(x, yPS7, yPS16, PS7_conc, PS16_conc, noPoints, nu, vol_solvent)
% from get_outputs in createFig4_Fig5.m / createFig6a.m / createFig6b.m
colours_MWDist = parula(1);
maxMW = max(x);
Nm = 1;
mass_polymer = PS7_conc;
weight_fraction_PS16 = 0;
y = yPS7*(1-weight_fraction_PS16)+yPS16*weight_fraction_PS16;
density_mass_PS = mass_polymer/vol_solvent;
[nPS7, MPS7, mol_polymer_densityPS7] = ...
    splitMWDist(false, colours_MWDist, x, y, noPoints, maxMW, Nm, nu, density_mass_PS);
mass_polymer = PS16_conc;
weight_fraction_PS16 = 1;
y = yPS7*(1-weight_fraction_PS16)+yPS16*weight_fraction_PS16;
density_mass_PS = mass_polymer/vol_solvent;
[nPS16, MPS16, mol_polymer_densityPS16] = ...
    splitMWDist(false, colours_MWDist, x, y, noPoints, maxMW, Nm, nu, density_mass_PS);
mol_polymer_density = mol_polymer_densityPS7 + mol_polymer_densityPS16;
mol_frac_PS7 = mol_polymer_densityPS7/mol_polymer_density;
n = [nPS16*(1-mol_frac_PS7), nPS7*mol_frac_PS7];
M = [MPS16, MPS7];
m = struct('n', n, 'M', M, 'molDensity', mol_polymer_density, 'nPS7', nPS7, 'MPS7', MPS7, ...
    'molDensityPS7', mol_polymer_densityPS7, 'nPS16', nPS16, 'MPS16', MPS16, ...
    'molDensityPS16', mol_polymer_densityPS16);
end

function pm = legacyModeProperties(n, M, mol_polymer_density, Ns, hstar, S, Ctau, ext, ...
    NA, kB, T, l, cinf, MW_monomer, K, alpha, eta_s)
% from get_outputs (createFig4_Fig5.m) and createFig2.m / createFig8.m
density = mol_polymer_density;
g = n'*density*kB*T*NA;
Ltrue = 2*l*M/MW_monomer;
R0 = sqrt(2*l^2.*M/MW_monomer*cinf);
Nk = Ltrue.^2./(2*cinf*l^2.*M/MW_monomer);
Nk = Nk';
if strcmp(ext, 'FENE-P')
    L = sqrt(3*((Ltrue'./R0').^2-1));
else
    L = Ltrue'./R0';
end
jvals = 1:Ns;
sigma = -1.40*hstar^0.78;
intrinsic_visc = K*M.^alpha;
lambda_m = 1/S*M.*intrinsic_visc*eta_s/(NA*kB*T)*Ctau;
lambda_ms = lambda_m'./jvals.^(2+sigma);
pm = struct('g', g, 'L', L, 'Nk', Nk, 'lambda_ms', lambda_ms, 'lambda_m', lambda_m, ...
    'Ltrue', Ltrue, 'R0', R0);
end

function dydt = pullODE(t, y, Nm, Ns, L, lambda,...
    Hfunc, Hderivfunc, zeta0, Lam0, V0, L0)
% odefun from Fig9_FastRetraction_CDD.m
    H = Hfunc(t*V0/L0);
    dH = Hderivfunc(t*V0/L0);
    z0 = ((zeta0-1)*zeta0)/(1-2*zeta0)^2;
    ep = (-3/(2*Lam0^2)*1./(H/(4*z0) - 1).*dH./H)*V0/L0;
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
