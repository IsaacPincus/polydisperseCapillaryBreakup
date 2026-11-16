function results = runTests(which)
%RUNTESTS Run the test suite.
%
%   runTests()           function and solver-physics tests (fast)
%   runTests('all')      also the unused-function and paper-result tests
%   runTests('paper')    paper-result tests only (slow: reruns every figure)
%   runTests('Fig3')     one paper figure (any name from paperFigures)
%
%   Generate the reference data first: generateReferenceData()

if nargin < 1
    which = 'functions';
end
testsDir = fileparts(mfilename('fullpath'));
addpath(fullfile(testsDir, 'helpers'));
setupTestPaths();

import matlab.unittest.TestSuite
switch which
    case 'functions'
        suite = [TestSuite.fromFile(fullfile(testsDir, 'FunctionTest.m')), ...
            TestSuite.fromFile(fullfile(testsDir, 'SolverPhysicsTest.m'))];
    case 'all'
        suite = [TestSuite.fromFile(fullfile(testsDir, 'FunctionTest.m')), ...
            TestSuite.fromFile(fullfile(testsDir, 'SolverPhysicsTest.m')), ...
            TestSuite.fromFile(fullfile(testsDir, 'UnusedFunctionTest.m')), ...
            TestSuite.fromFile(fullfile(testsDir, 'PaperResultTest.m'))];
    case 'paper'
        suite = TestSuite.fromFile(fullfile(testsDir, 'PaperResultTest.m'));
    otherwise
        suite = TestSuite.fromFile(fullfile(testsDir, 'PaperResultTest.m'), ...
            'ParameterProperty', 'figureName', 'ParameterName', which);
end
runner = matlab.unittest.TestRunner.withTextOutput('Verbosity', 2);
results = runner.run(suite);
disp(table(results))
end
