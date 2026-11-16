function r = computeFig4_Fig5()
%COMPUTEFIG4_FIG5 Single-species thinning: simulation vs experiment (Figs 4, 5).
%
%   Fig 4a: PS16 at 5-200 wppm; Fig 4b: PS7 at 200-1000 wppm; Fig 5:
%   tau_EC/tau_Z against c/c*, and (inset) f_min tau_EC/tau_Z with f_min
%   the FENE-PM factor of that species at the strain-rate minimum.
%   Constant drag, Ns = 10.

p = defaultParameters();
Ns = 10;
cstarPS7 = 3.92e3;    % wppm
cstarPS16 = 2.59e3;   % wppm
mw = loadMWDistributions(p);

% --- Fig 4a: PS16 only --------------------------------------------------
s.cPS16 = [5, 10, 20, 50, 100, 200];
s.cPS7 = zeros(size(s.cPS16));
[s.t_exp, s.D_exp] = readThinningCurves('PS16', numel(s.cPS16));
s.sim = simulateConcentrationSweep(p, mw, s.cPS7, s.cPS16, Ns, 'constant');
s = analyseExperiments(s, num2cell(s.t_exp, 1), num2cell(s.D_exp, 1), ...
    'Threshold', 5, 'FitRadius', true, 'StartPoint', true);
s.tau_EC_ratio = s.sim.tau_EC/s.sim.lambda_m(1);
s.tau_EC_exp_ratio = s.tau_EC_exp/s.sim.lambda_m(1);
s.f_min_ratio = s.sim.f_min(:,1)'.*s.tau_EC_ratio;
s.M_implied = (s.tau_EC_exp/(p.K*p.eta_s)*(p.NA*p.kB*p.T*p.S)).^(2/3);
r.PS16 = s;

% --- Fig 4b: PS7 only ---------------------------------------------------
s = struct();
s.cPS7 = [200, 500, 1000];
s.cPS16 = zeros(size(s.cPS7));
s.t_exp = cell(1, 3);
s.D_exp = cell(1, 3);
for ii = 1:3
    sheet = sprintf('%dppmPS7+PS16', s.cPS7(ii));
    [s.t_exp{ii}, s.D_exp{ii}] = readThinningCurves(sheet, 1, false);
end
s.sim = simulateConcentrationSweep(p, mw, s.cPS7, s.cPS16, Ns, 'constant');
s = analyseExperiments(s, s.t_exp, s.D_exp, 'Threshold', 3, 'FitRadius', true);
s.tau_EC_ratio = s.sim.tau_EC/s.sim.lambda_m(2);
s.tau_EC_exp_ratio = s.tau_EC_exp/s.sim.lambda_m(2);
s.f_min_ratio = s.sim.f_min(:,2)'.*s.tau_EC_ratio;
r.PS7 = s;

% --- Fig 5: c/c* ----------------------------------------------------------
r.cr7 = r.PS7.cPS7/cstarPS7;
r.cr16 = r.PS16.cPS16/cstarPS16;
end

function s = analyseExperiments(s, t, D, varargin)
for ii = 1:numel(t)
    e = experimentalRelaxationTime(t{ii}, D{ii}, varargin{:});
    s.tau_ep_min_exp(ii) = e.tau_ep_min;
    s.tau_EC_exp(ii) = e.tau_EC;
    s.tavg_exp{ii} = e.tavg;
    s.ep_exp{ii} = e.ep;
end
end
