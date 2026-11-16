function d = referenceDirectory()
%REFERENCEDIRECTORY Folder holding the reference .mat files.
d = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'reference');
end
