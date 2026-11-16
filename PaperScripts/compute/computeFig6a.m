function r = computeFig6a()
%COMPUTEFIG6A Blend thinning curves with conformation-dependent drag (Fig 6a).
%
%   PS16 at 0-200 wppm added to 0, 200 and 500 wppm PS7, Ns = 10.
%   The paper panel is the 200 wppm PS7 set.

p = defaultParameters();
Ns = 10;
mw = loadMWDistributions(p);

% cPS7, spreadsheet sheet, cPS16 values, experimental fit options
sets = {
    0,   'PS16',           [5,10,20,50,100,200],   {'Threshold', 5, 'StartPoint', true}
    200, '200ppmPS7+PS16', [0,5,10,20,50,100,200], {'Threshold', 3}
    500, '500ppmPS7+PS16', [0,5,10,20,50,100,200], {'Threshold', 5}
    };
r.blocks = cell(1, size(sets, 1));
for k = 1:size(sets, 1)
    r.blocks{k} = blendSet(p, mw, Ns, sets{k,:});
end
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
    b.tavg_exp{ii} = e.tavg;
    b.ep_exp{ii} = e.ep;
end
end
