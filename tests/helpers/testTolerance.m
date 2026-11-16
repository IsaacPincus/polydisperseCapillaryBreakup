function relTol = testTolerance()
%TESTTOLERANCE Relative tolerance used when comparing with the original code.
%   The refactor keeps every arithmetic operation in its original order, so
%   on the same MATLAB installation results are normally bit-identical; the
%   tests report this. The tolerance only absorbs round-off differences.
relTol = 1e-9;
end
