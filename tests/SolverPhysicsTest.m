classdef SolverPhysicsTest < matlab.unittest.TestCase
    %SOLVERPHYSICSTEST solveFilamentThinning satisfies the paper's equations.
    %
    %   Covers what the comparison with the original code cannot: multi-
    %   mode chains (Eq. 22 with 3Ns) and conformation-dependent drag.
    %   See helpers/checkSolverPhysics.

    methods (TestClassSetup)
        function paths(testCase)
            addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
            setupTestPaths();
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
                'MATLAB:nearlySingularMatrix'));
        end
    end

    methods (Test, TestTags = {'Functions'})
        function satisfiesPaperEquations(testCase)
            problems = checkSolverPhysics();
            testCase.verifyEmpty(problems, strjoin(problems, newline));
        end
    end
end
