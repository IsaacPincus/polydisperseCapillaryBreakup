function r = computeFig6b()
%COMPUTEFIG6B tau_EC of PS7 + PS16 blends vs PS16 concentration (Fig 6b).
%
%   PS16 at 0-200 wppm in 0, 200, 500 and 1000 wppm PS7, conformation-
%   dependent drag, Ns = 30. Also loads the Haward et al. data shown in
%   the figure and the tau_EC/tau_Z ratios plotted against c/c*.

p = defaultParameters();
Ns = 30;
cstarPS7 = 3.92e3;    % wppm
cstarPS16 = 2.59e3;   % wppm
mw = loadMWDistributions(p);

% cPS7, spreadsheet sheet, cPS16 values, experimental fit options
sets = {
    0,    'PS16',            [5,10,20,50,100,200],   {'StartPoint', true}
    200,  '200ppmPS7+PS16',  [0,5,10,20,50,100,200], {'PositiveRatesOnly', true}
    500,  '500ppmPS7+PS16',  [0,5,10,20,50,100,200], {}
    1000, '1000ppmPS7+PS16', [0,5,10,20,50,100,200], {}
    };
r.blocks = cell(1, size(sets, 1));
for k = 1:size(sets, 1)
    r.blocks{k} = blendSet(p, mw, Ns, sets{k,:});
end

% tau_EC/tau_Z for the single-species cases (Fig 5-style comparison)
b = r.blocks{1};
r.tau_EC_ratio_PS16 = b.sim.tau_EC/b.sim.lambda_m(1);
r.tau_EC_exp_ratio_PS16 = b.tau_EC_exp/b.sim.lambda_m(1);
for k = 2:4
    b = r.blocks{k};
    r.tau_EC_ratio_PS7(k-1) = b.sim.tau_EC(1)/b.sim.lambda_m(2);
    r.tau_EC_exp_ratio_PS7(k-1) = b.tau_EC_exp(1)/b.sim.lambda_m(2);
end
r.cr7 = [200, 500, 1000]/cstarPS7;
r.cr16 = [5, 10, 20, 50, 100, 200]/cstarPS16;

% Haward et al.: columns (c_PS16, tau_EC [ms]) for 0/200/500/1000 wppm PS7
r.haward = readDataFile('Fig5_Hawardetal.csv');
end

function b = blendSet(p, mw, Ns, cPS7, sheet, cPS16, fitOptions)
b.cPS16 = cPS16;
b.cPS7 = cPS7*ones(size(cPS16));
[b.t_exp, b.D_exp] = readThinningCurves(sheet, numel(cPS16));
b.sim = simulateConcentrationSweep(p, mw, b.cPS7, b.cPS16, Ns, 'conformationDependent');
for ii = 1:numel(cPS16)
    e = experimentalRelaxationTime(b.t_exp(:,ii), b.D_exp(:,ii), fitOptions{:});
    b.tau_ep_min_exp(ii) = e.tau_ep_min;
    b.tau_EC_exp(ii) = e.tau_EC;
end
end
