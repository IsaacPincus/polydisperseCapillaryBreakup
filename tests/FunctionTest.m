classdef FunctionTest < matlab.unittest.TestCase
    %FUNCTIONTEST Functions reproduce the original code where they should.
    %
    %   Each test runs a function from Functions/ on the inputs stored in
    %   tests/reference/functions.mat and compares the outputs with what
    %   the original code produced on the same inputs (see
    %   generateReferenceData and helpers/legacyFunctionReferences). Cases
    %   where the paper-consistent code deliberately differs from the
    %   original are covered by SolverPhysicsTest instead.
    %
    %   Run with runTests or: runtests('FunctionTest')

    properties (TestParameter)
        caseName = {'inertialDensity', 'ubcal', 'zetacal', 'loadMWDistributions', ...
            'splitMWDist', 'twoSpeciesModes', 'blendModes', 'modeProperties', ...
            'linearViscoelasticModuli', 'packConformation', 'readThinningCurves', ...
            'solverConstantDrag', 'solverConstantDragSingleMode', ...
            'solverInitialState', 'simulatedRelaxationTime', ...
            'experimentalRelaxationTime', 'plateRetraction'};
    end

    properties
        cases
    end

    methods (TestClassSetup)
        function loadReferences(testCase)
            addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
            setupTestPaths();
            file = fullfile(referenceDirectory(), 'functions.mat');
            testCase.assertTrue(isfile(file), ...
                'No reference data: run generateReferenceData() first.');
            data = load(file);
            testCase.cases = data.cases;
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
                'MATLAB:nearlySingularMatrix'));
        end
    end

    methods (Test, TestTags = {'Functions'})
        function matchesOriginal(testCase, caseName)
            testCase.assertTrue(isfield(testCase.cases, caseName), ...
                sprintf('No reference for %s: rerun generateReferenceData().', caseName));
            c = testCase.cases.(caseName);
            actual = currentFunctionOutputs(caseName, c.inputs);
            [problems, identical] = compareValues(actual, c.outputs, caseName, testTolerance());
            testCase.verifyEmpty(problems, strjoin(problems, newline));
            if isempty(problems) && ~identical
                testCase.log(1, sprintf('%s agrees within tolerance but is not bit-identical.', caseName));
            end
        end
    end
end
