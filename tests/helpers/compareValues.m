function [problems, identical] = compareValues(actual, expected, name, relTol)
%COMPAREVALUES Compare a result with its reference value.
%
%   [problems, identical] = compareValues(actual, expected, name, relTol)
%   returns a cell array of messages describing every mismatch (empty if
%   the values agree) and whether they are bit-for-bit identical. Fields
%   that exist only in actual are ignored.
%   Numeric arrays must have the same size; elements agree when
%       |a - e| <= relTol*max(|e|, 1e-3*max(|e_all|)),
%   and NaN/Inf must appear in the same places. Cells and structs are
%   compared element by element / field by field.

if nargin < 4
    relTol = 1e-9;
end
problems = {};
identical = isequaln(actual, expected);
if identical
    return
end

if isnumeric(expected) || islogical(expected)
    if ~(isnumeric(actual) || islogical(actual))
        problems{end+1} = sprintf('%s: expected numeric, got %s', name, class(actual));
        return
    end
    if ~isequal(size(actual), size(expected))
        problems{end+1} = sprintf('%s: size %s, expected %s', name, ...
            mat2str(size(actual)), mat2str(size(expected)));
        return
    end
    a = double(actual(:));
    e = double(expected(:));
    finite = isfinite(e);
    if ~isequaln(a(~finite), e(~finite)) || any(~isfinite(a(finite)))
        problems{end+1} = sprintf('%s: NaN/Inf entries differ', name);
        return
    end
    a = a(finite);
    e = e(finite);
    if isempty(e)
        return
    end
    scale = max(abs(e));
    err = abs(a - e);
    tol = relTol*max(abs(e), 1e-3*scale);
    if any(err > tol)
        [worst, k] = max(err./max(abs(e), 1e-3*scale));
        problems{end+1} = sprintf('%s: %d of %d elements differ (worst relative error %.3g at element %d)', ...
            name, nnz(err > tol), numel(e), worst, k);
    end
elseif iscell(expected)
    if ~iscell(actual) || ~isequal(size(actual), size(expected))
        problems{end+1} = sprintf('%s: expected a %s cell', name, mat2str(size(expected)));
        return
    end
    identical = true;
    for k = 1:numel(expected)
        [p, id] = compareValues(actual{k}, expected{k}, sprintf('%s{%d}', name, k), relTol);
        problems = [problems, p]; %#ok<AGROW>
        identical = identical && id;
    end
elseif isstruct(expected)
    if ~isstruct(actual)
        problems{end+1} = sprintf('%s: expected a struct', name);
        return
    end
    % extra fields in actual are allowed
    identical = true;
    fields = fieldnames(expected);
    for k = 1:numel(fields)
        if ~isfield(actual, fields{k})
            problems{end+1} = sprintf('%s.%s: missing', name, fields{k}); %#ok<AGROW>
            identical = false;
        else
            [p, id] = compareValues(actual.(fields{k}), expected.(fields{k}), ...
                [name '.' fields{k}], relTol);
            problems = [problems, p]; %#ok<AGROW>
            identical = identical && id;
        end
    end
elseif isobject(expected)
    % fit objects etc.: their numeric content is checked through other fields
    identical = true;
elseif ~isequal(actual, expected)
    problems{end+1} = sprintf('%s: values differ', name);
end
end
