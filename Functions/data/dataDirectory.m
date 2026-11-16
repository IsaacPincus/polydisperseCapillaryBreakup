function d = dataDirectory()
%DATADIRECTORY Absolute path of the project's Data folder.

projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
d = fullfile(projectRoot, 'Data');
end
