function v = packConformation(Az, Ar)
%PACKCONFORMATION Interleave Nm-by-Ns Az and Ar into a column vector.
%
%   v = [Az(1,1) Ar(1,1) Az(1,2) Ar(1,2) ... Az(Nm,Ns) Ar(Nm,Ns)]'.
%   Inverse of unpackConformation.

v = reshape(permute(cat(3, Az, Ar), [3 2 1]), [], 1);
end
