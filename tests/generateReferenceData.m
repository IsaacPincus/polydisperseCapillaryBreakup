function generateReferenceData(varargin)
%GENERATEREFERENCEDATA Record the reference results the tests compare against.
%
%   generateReferenceData()   both kinds of reference (run once)
%
%   Functions: the ORIGINAL functions (git commit LEGACY_COMMIT, extracted
%   with 'git archive', or a folder given with 'LegacyDir') are run on fixed
%   inputs -> tests/reference/functions.mat, unusedFunctions.mat. The
%   cases are chosen where the paper-consistent code should agree with the
%   original; see helpers/legacyFunctionReferences.
%
%   Figures: the CURRENT computeFig* functions are run and their results
%   saved -> tests/reference/<figure>.mat. PaperResultTest then guards
%   these numbers against unintended changes. Regenerate a figure after a
%   deliberate change to its physics or parameters.
%
%   Options
%     'Figures'   names from paperFigures to (re)generate, default all,
%                 e.g. {'Fig3', 'Fig7'}; {} skips the figures
%     'Functions' true (default) to (re)generate the function references
%     'LegacyDir' folder holding the original Functions/ and Data/
%     'Commit'    git commit with the original code (default LEGACY_COMMIT)
%
%   Running every figure takes as long as running them by hand (Fig 8
%   dominates).

LEGACY_COMMIT = '24f47bb';

testsDir = fileparts(mfilename('fullpath'));
root = fileparts(testsDir);
addpath(fullfile(testsDir, 'helpers'));
figs = paperFigures();

ip = inputParser;
ip.addParameter('Figures', {figs.name});
ip.addParameter('Functions', true);
ip.addParameter('LegacyDir', '');
ip.addParameter('Commit', LEGACY_COMMIT);
ip.parse(varargin{:});
opt = ip.Results;
if ischar(opt.Figures)
    opt.Figures = {opt.Figures};
end

refDir = fullfile(testsDir, 'reference');
if ~exist(refDir, 'dir')
    mkdir(refDir);
end
warning('off', 'MATLAB:nearlySingularMatrix');

if opt.Functions
    legacyDir = opt.LegacyDir;
    if isempty(legacyDir)
        legacyDir = extractLegacyCode(root, opt.Commit);
        removeLegacy = onCleanup(@() rmdir(legacyDir, 's'));
    end
    fprintf('Function references from the original code in %s ...\n', legacyDir);
    meta = struct('source', 'original', 'commit', opt.Commit, ...
        'created', datestr(now), 'version', version()); %#ok<TNOW1,DATST>
    % only the original code may be on the path: it shares function names
    oldPath = path();
    restorePath = onCleanup(@() path(oldPath));
    rmpath(genpath(fullfile(root, 'Functions')));
    rmpath(genpath(fullfile(root, 'PaperScripts')));
    addpath(fullfile(legacyDir, 'Functions'), fullfile(legacyDir, 'Data'));
    [cases, unusedCases] = legacyFunctionReferences(legacyDir); %#ok<ASGLU>
    save(fullfile(refDir, 'functions.mat'), 'cases', 'meta');
    save(fullfile(refDir, 'unusedFunctions.mat'), 'unusedCases', 'meta');
    clear restorePath   % restore the path before deleting the extracted code
    if exist('removeLegacy', 'var')
        clear removeLegacy
    end
end

setupTestPaths();
oldVisible = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');
restoreVisible = onCleanup(@() set(0, 'DefaultFigureVisible', oldVisible)); %#ok<NASGU>
for f = 1:numel(figs)
    if ~any(strcmp(figs(f).name, opt.Figures))
        continue
    end
    fprintf('%s ...\n', figs(f).name);
    tic
    results = feval(figs(f).compute); %#ok<NASGU>
    meta = struct('source', 'current', 'created', datestr(now), 'version', version()); %#ok<TNOW1,DATST,NASGU>
    save(fullfile(refDir, [figs(f).name '.mat']), 'results', 'meta');
    fprintf('  %.0f s\n', toc);
end
fprintf('References saved in %s\n', refDir);
end

function legacyDir = extractLegacyCode(root, commit)
legacyDir = tempname;
mkdir(legacyDir);
zipFile = [legacyDir '.zip'];
command = sprintf('git -C "%s" archive --format=zip -o "%s" %s Functions Data', ...
    root, zipFile, commit);
[status, output] = system(command);
if status ~= 0
    error('generateReferenceData:git', ['Could not extract commit %s with git:\n%s\n' ...
        'Check out the original code into a folder and pass it with ''LegacyDir''.'], ...
        commit, output);
end
unzip(zipFile, legacyDir);
delete(zipFile);
end
