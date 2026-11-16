function r = simulatedRelaxationTime(t, a, ep, varargin)
%SIMULATEDRELAXATIONTIME Elasto-capillary relaxation time of a simulation.
%
%   r = simulatedRelaxationTime(t, a, ep, Name, Value)
%
%   Implements Fig. 3 / Eq. 32 of the paper. The critical time t_c is the
%   first local maximum of ep (or 0). ep_min is the first local minimum of
%   ep after t_c (or the global minimum). The EC regime is every point with
%   t > t_c and ep < Threshold*ep_min (within 10% of ep_min by default),
%   where a = a1 exp(-(t - t_c)/(3 tau_EC)) is fitted.
%
%   Options (defaults in brackets)
%     'Threshold'    [1.1]   EC regime: ep < Threshold*ep_min
%     'ShiftTime'    [true]  fit against t - t_c (false: against t;
%                            same tau_EC up to round-off)
%     'CriticalTime' ['firstMaximum'] or 'zero'
%     'Minimum'      ['localAfterTc'] or 'global'
%
%   Output fields: t_c, ep_min, ep_min_index, t_ep_min (time of ep_min),
%   tau_ep_min (2/(3 ep_min)), indices_EC, threshold, fitobject, tau_EC.

opt = parseOptions(struct('Threshold', 1.1, 'ShiftTime', true, ...
    'CriticalTime', 'firstMaximum', 'Minimum', 'localAfterTc'), varargin{:});

switch opt.CriticalTime
    case 'firstMaximum'
        t_max_ep_all = t(islocalmax(ep)==1);
        r.t_c = t_max_ep_all(1);
    case 'zero'
        r.t_c = 0;
end

switch opt.Minimum
    case 'localAfterTc'
        r.ep_min_index = find((islocalmin(ep)==1)&(t>r.t_c), 1);
        r.ep_min = ep(r.ep_min_index);
    case 'global'
        [r.ep_min, r.ep_min_index] = min(ep);
end
r.t_ep_min = t(r.ep_min_index);
r.tau_ep_min = 2/(3*r.ep_min);

r.threshold = opt.Threshold;
r.indices_EC = find((ep<opt.Threshold*r.ep_min)&(t>r.t_c));

tFit = t(r.indices_EC);
if opt.ShiftTime
    tFit = tFit - r.t_c;
end
[xxx, yyy] = prepareCurveData(tFit, a(r.indices_EC));
r.fitobject = fit(xxx, yyy, 'exp1');
r.tau_EC = -1./(3*r.fitobject.b);
end
