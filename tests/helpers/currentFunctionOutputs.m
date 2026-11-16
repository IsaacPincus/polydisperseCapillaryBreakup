function out = currentFunctionOutputs(caseName, in)
%CURRENTFUNCTIONOUTPUTS Run the REFACTORED code on a function-test case.
%
%   out = currentFunctionOutputs(caseName, inputs) returns outputs laid
%   out like the references from legacyFunctionReferences, so the two can
%   be compared with compareValues.

switch caseName
    case 'inertialDensity'
        out.rhoEff = inertialDensity(in.rho, in.alpha_C1);

    case 'ubcal'
        out.Ub = cellfun(@(h) ubcal(h, in.nkperb), in.hstar, 'UniformOutput', false);

    case 'zetacal'
        out = struct('zet', nan(size(in.Q)), 'G', nan(size(in.Q)), ...
            'zb', nan(size(in.Q)), 'Nb', nan(size(in.Q)));
        for k = 1:numel(in.Q)
            [out.zet(k), G, zb, Nb] = zetacal(in.Q(k), 1, in.phieq(k), in.nk, in.hsk(k), in.himodel(k));
            out.G(k) = G; out.zb(k) = zb; out.Nb(k) = Nb;
        end

    case 'loadMWDistributions'
        mw = loadMWDistributions(defaultParameters());
        out = struct('x', mw.x, 'yPS7', mw.yPS7, 'yPS16', mw.yPS16);

    case 'splitMWDist'
        out.results = cell(1, numel(in.N));
        for k = 1:numel(in.N)
            [n, M, mol] = splitMWDist(in.x, in.y, in.N{k}, 0.5, in.density);
            out.results{k} = struct('n', n, 'M', M, 'molDensity', mol);
        end

    case 'twoSpeciesModes'
        p = defaultParameters();
        mw = loadMWDistributions(p);
        out.results = cell(1, size(in.conc, 1));
        for k = 1:size(in.conc, 1)
            out.results{k} = twoSpeciesModes(mw, in.conc(k,1), in.conc(k,2), p);
        end

    case 'blendModes'
        p = defaultParameters();
        mw = loadMWDistributions(p);
        out.results = cell(1, size(in.conc, 1));
        for k = 1:size(in.conc, 1)
            out.results{k} = blendModes(mw, in.conc(k,1), in.conc(k,2), in.conc(k,3), p);
        end

    case 'modeProperties'
        p = defaultParameters();
        out.results = cell(1, size(in.variants, 1));
        for k = 1:size(in.variants, 1)
            [Ns, hstar, S] = in.variants{k,:};
            out.results{k} = modeProperties(in.modes, p, Ns, 'hstar', hstar, 'S', S);
        end

    case 'linearViscoelasticModuli'
        [out.Gp, out.Gpp] = linearViscoelasticModuli(in.pm, in.freq);

    case 'packConformation'
        [out.Az, out.Ar] = unpackConformation(in.v, in.Nm, in.Ns);
        assert(isequal(packConformation(out.Az, out.Ar), in.v), ...
            'packConformation does not invert unpackConformation');

    case {'solverConstantDrag', 'solverConstantDragSingleMode'}
        out = solveFilamentThinning(in.p, in.modes);

    case 'solverInitialState'
        out = solveFilamentThinning(in.p, in.modes, ...
            'InitialConformation', in.conformation, 'InitialStrainRate', in.epInitial);

    case 'simulatedRelaxationTime'
        % the original scripts' windows (the paper's default is 10%)
        out.getOutputs = simulatedRelaxationTime(in.t, in.a, in.ep, ...
            'Threshold', 1.2, 'ShiftTime', false);
        r = simulatedRelaxationTime(in.t, in.a, in.ep, 'Threshold', 1.2);
        out.shifted = struct('t_c', r.t_c, 't_ep_min', r.t_ep_min, 'tau_EC', r.tau_EC, ...
            'fitA', r.fitobject.a);
        out.threshold2 = simulatedRelaxationTime(in.t, in.a, in.ep, ...
            'Threshold', 2, 'ShiftTime', false);
        out.globalMinimum = simulatedRelaxationTime(in.t, in.a, in.ep, ...
            'Threshold', 1.2, 'CriticalTime', 'zero', 'Minimum', 'global');

    case 'readThinningCurves'
        [out.t16, out.D16] = readThinningCurves('PS16', 6);
        [out.t7, out.D7] = readThinningCurves('200ppmPS7+PS16', 1, false);
        [out.tMix, out.DMix] = readThinningCurves('200ppmPS7+PS16', 7);

    case 'experimentalRelaxationTime'
        out.radiusStartPoint = experimentalRelaxationTime(in.t16, in.D16, ...
            'Threshold', 5, 'FitRadius', true, 'StartPoint', true);
        out.radiusThreshold3 = experimentalRelaxationTime(in.t7, in.D7, ...
            'Threshold', 3, 'FitRadius', true);
        out.diameterPositive = experimentalRelaxationTime(in.tMix, in.DMix, ...
            'PositiveRatesOnly', true);

    case 'plateRetraction'
        out = plateRetraction(in.pm, in.V0, in.geom);

    otherwise
        error('currentFunctionOutputs:case', 'No adapter for case "%s".', caseName);
end
end
