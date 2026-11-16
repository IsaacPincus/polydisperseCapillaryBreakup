function opt = parseOptions(defaults, varargin)
%PARSEOPTIONS Merge Name,Value pairs into a struct of defaults.
%
%   opt = parseOptions(defaults, 'Name1', value1, ...) errors on names
%   that are not fields of defaults (case-sensitive).

opt = defaults;
if mod(numel(varargin), 2) ~= 0
    error('parseOptions:pairs', 'Options must be given as Name,Value pairs.');
end
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~isfield(defaults, name)
        error('parseOptions:unknown', 'Unknown option "%s".', name);
    end
    opt.(name) = varargin{k+1};
end
end
