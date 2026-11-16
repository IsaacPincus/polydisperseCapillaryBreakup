classdef PaperResultTest < matlab.unittest.TestCase
    %PAPERRESULTTEST Paper results have not changed since they were recorded.
    %
    %   Reruns each computeFig* function and compares every number in its
    %   results with tests/reference/<figure>.mat, recorded from the same
    %   code by generateReferenceData. Any difference beyond round-off means
    %   a change altered the paper's results. These tests are slow: they
    %   rerun every simulation (Fig 8 takes longest). See runTests for
    %   running one figure.

    properties (TestParameter)
        figureName = struct('Fig1', 'Fig1', 'Fig2', 'Fig2', 'Fig3', 'Fig3', ...
            'Fig4_Fig5', 'Fig4_Fig5', 'Fig6a', 'Fig6a', 'Fig6b', 'Fig6b', ...
            'Fig7', 'Fig7', 'Fig8', 'Fig8', 'Fig9_Fig10', 'Fig9_Fig10');
    end

    methods (TestClassSetup)
        function paths(testCase)
            addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
            setupTestPaths();
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
                'MATLAB:nearlySingularMatrix'));
        end
    end

    methods (Test, TestTags = {'Paper'})
        function matchesRecordedResults(testCase, figureName)
            file = fullfile(referenceDirectory(), [figureName '.mat']);
            testCase.assumeTrue(isfile(file), sprintf( ...
                'No reference for %s: run generateReferenceData(''Figures'', {''%s''}).', ...
                figureName, figureName));
            data = load(file);
            expected = data.results;
            actual = feval(['compute' figureName]);
            [problems, identical] = compareValues(actual, expected, figureName, testTolerance());
            testCase.verifyEmpty(problems, strjoin(problems, newline));
            if isempty(problems) && ~identical
                testCase.log(1, sprintf('%s agrees within tolerance but is not bit-identical.', figureName));
            end
        end
    end
end
