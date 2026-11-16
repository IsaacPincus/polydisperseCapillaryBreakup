classdef UnusedFunctionTest < matlab.unittest.TestCase
    %UNUSEDFUNCTIONTEST Functions in Functions/unused still behave as before.
    %
    %   No paper script calls these functions. They were moved, unchanged,
    %   to Functions/unused; this test guards them until it is decided
    %   whether to keep them. Delete this file together with that folder.

    properties (TestParameter)
        caseName = {'getInitialStretch', 'getTimes', 'fitGaussiansToMWDist', ...
            'getStretchExperiments', 'getStretchExperimentsStressBalanced', ...
            'filamentThinningFENEPM_inertia_initialEXP', ...
            'filamentThinningFENEPM_inertia_simple_CDD'};
    end

    properties
        cases
    end

    methods (TestClassSetup)
        function loadReferences(testCase)
            addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
            setupTestPaths();
            file = fullfile(referenceDirectory(), 'unusedFunctions.mat');
            testCase.assertTrue(isfile(file), ...
                'No reference data: run generateReferenceData() first.');
            data = load(file);
            testCase.cases = data.unusedCases;
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
                'MATLAB:nearlySingularMatrix'));
        end
    end

    methods (Test, TestTags = {'Unused'})
        function matchesOriginal(testCase, caseName)
            c = testCase.cases.(caseName);
            testCase.assumeEmpty(c.error, ['The original function failed on the test input: ' c.error]);
            in = c.inputs;
            names = fieldnames(c.outputs)';
            out = cell(1, numel(names));
            switch caseName
                case 'getInitialStretch'
                    [out{:}] = getInitialStretch(in.L0, in.L1, in.R0);
                case 'getTimes'
                    s = in.sol;
                    [out{:}] = getTimes(s.a, s.Az, s.Ar, s.f, s.stress, s.ep, s.Wi, s.t, ...
                        s.te, s.ye, s.ie, in.chi, in.mu, in.N, in.L, in.lambda, in.g, in.a0, false);
                case 'fitGaussiansToMWDist'
                    [out{:}] = fitGaussiansToMWDist(false, parula(2), in.x, in.y, ...
                        numel(in.x), max(in.x), 2, 0.5, in.density);
                case 'getStretchExperiments'
                    [out{:}] = getStretchExperiments(in.texp, in.Dexp, 2, 3, in.L, in.lambda);
                case 'getStretchExperimentsStressBalanced'
                    [out{:}] = getStretchExperimentsStressBalanced(in.p.chi, in.p.X, ...
                        in.p.eta_s, in.p.rho, 2, 3, in.modes.L, in.modes.lambda_ms, ...
                        in.modes.g, in.texp, in.Dexp);
                case 'filamentThinningFENEPM_inertia_initialEXP'
                    [out{:}] = filamentThinningFENEPM_inertia_initialEXP(in.p.chi, in.p.X, ...
                        in.p.eta_s, in.p.rho, 2, 3, in.modes.L, in.modes.lambda_ms, ...
                        in.modes.g, in.p.cutoff, in.p.tmax, in.texp, in.Dexp);
                case 'filamentThinningFENEPM_inertia_simple_CDD'
                    [out{:}] = filamentThinningFENEPM_inertia_simple_CDD(in.p.chi, in.p.X, ...
                        in.p.eta_s, in.p.rho, 2, 3, in.modes.L, in.modes.Nk, in.zet_rat, ...
                        in.modes.lambda_ms, in.modes.g, in.p.a0, in.p.cutoff, in.p.tmax);
            end
            actual = cell2struct(out, names, 2);
            problems = compareValues(actual, c.outputs, caseName, testTolerance());
            testCase.verifyEmpty(problems, strjoin(problems, newline));
        end
    end
end
