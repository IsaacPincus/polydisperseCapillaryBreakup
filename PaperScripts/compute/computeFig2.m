function r = computeFig2()
%COMPUTEFIG2 Linear viscoelastic moduli of PS7, PS16 and their blend (Fig 2).
%
%   Each species is one MW mode with 50 Zimm sub-modes. Also returns the
%   Ohnesorge and elasto-capillary numbers and viscosity ratio of Table I.

p = defaultParameters();
Ns = 50;
mw = loadMWDistributions(p);

r.freq = logspace(-1, 3, 100);   % rad/s

r.PS7 = species(mw, p, 1000, 0, Ns, r.freq);
r.PS7.exp = readSAOS('PS7_1000ppm_SAOS.xlsx', p.eta_s);
r.PS16 = species(mw, p, 0, 200, Ns, r.freq);
r.PS16.exp = readSAOS('PS16_200ppm_SAOS.xlsx', p.eta_s);

r.mix.Gp = r.PS7.Gp + r.PS16.Gp;
r.mix.Gpp = r.PS7.Gpp + r.PS16.Gpp;
SAOS = readDataFile('SAOS_mix_1000PS7_200PS16.csv');
r.mix.exp.freq = SAOS(:,1);
r.mix.exp.Gp = SAOS(:,2);
r.mix.exp.Gpp = SAOS(:,3) - p.eta_s*r.mix.exp.freq;
end

function s = species(mw, p, cPS7, cPS16, Ns, freq)
modes = blendModes(mw, cPS7, cPS16, 1, p);
pm = modeProperties(modes, p, Ns);
s.n = modes.n;
s.M = modes.M;
s.g = pm.g;
s.L = pm.L;
s.lambda_ms = pm.lambda_ms;

Rg = pm.R0/sqrt(6);
s.cstar = modes.M/(p.NA*(2*Rg)^3);
eta_l = p.eta_s + sum(pm.lambda_ms.*pm.g);
s.Oh = eta_l/sqrt(p.rhoFluid*p.chi*p.a0);   % Ohnesorge number (Eq. 12a)
s.Ec = p.chi/(pm.g*p.a0);            % elasto-capillary number (Eq. 12b)
s.beta = p.eta_s/eta_l;
[s.Gp, s.Gpp] = linearViscoelasticModuli(pm, freq);
end

function e = readSAOS(fileName, eta_s)
% SAOS data with G'' corrected for the solvent contribution
SAOS = readDataFile(fileName);
e.freq = SAOS(:,1);
e.Gp = SAOS(:,2);
e.Gp_std = SAOS(:,3);
e.Gpp = SAOS(:,4) - eta_s*e.freq;
e.Gpp_std = SAOS(:,5);
end
