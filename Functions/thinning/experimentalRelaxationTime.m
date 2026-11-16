function r = experimentalRelaxationTime(t, D, varargin)
%EXPERIMENTALRELAXATIONTIME Relaxation time from a measured CaBER curve.
%
%   r = experimentalRelaxationTime(t, D, Name, Value) takes t - t_c [s]
%   and diameter D [mm] of one experiment, computes the strain rate
%   ep = -2/y dy/dt at interval midpoints (y = D or D/2), and fits
%   y = A exp(b t) over t > 0 with ep < Threshold*min(ep): tau_EC = -1/(3b).
%
%   Options (defaults in brackets)
%     'Threshold'         [5]
%     'FitRadius'         [false]  fit a = D/2 instead of D
%     'StartPoint'        [false]  seed the fit at [0.1, -1/(3 tau_ep_min)]
%     'PositiveRatesOnly' [false]  take min(ep) over ep > 0 only
%
%   Output fields: tavg, yavg, ep, ep_min, tau_ep_min, indices_EC,
%   fitobject, tau_EC.

opt = parseOptions(struct('Threshold', 5, 'FitRadius', false, ...
    'StartPoint', false, 'PositiveRatesOnly', false), varargin{:});

if opt.FitRadius
    y = D/2;
else
    y = D;
end
r.tavg = (t(2:end) + t(1:end-1))/2;
r.yavg = (y(2:end) + y(1:end-1))/2;
r.ep = -2./r.yavg.*diff(y)./diff(t);

window = r.tavg>0;
if opt.PositiveRatesOnly
    window = window & r.ep>0;
end
r.ep_min = min(r.ep(window));
r.tau_ep_min = 2/(3*r.ep_min);

r.indices_EC = find((r.ep<opt.Threshold*r.ep_min)&(r.tavg>0));
[xxx, yyy] = prepareCurveData(r.tavg(r.indices_EC), r.yavg(r.indices_EC));
if opt.StartPoint
    r.fitobject = fit(xxx, yyy, 'exp1', 'StartPoint', [0.1, -1/(3*r.tau_ep_min)]);
else
    r.fitobject = fit(xxx, yyy, 'exp1');
end
r.tau_EC = -1./(3*r.fitobject.b);
end
