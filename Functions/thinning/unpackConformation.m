function [Az, Ar] = unpackConformation(v, Nm, Ns)
%UNPACKCONFORMATION Split a state vector into Nm-by-Ns Az and Ar.
%
%   The vector is ordered [Az(1,1) Ar(1,1) Az(1,2) Ar(1,2) ... Az(Nm,Ns) Ar(Nm,Ns)]
%   (sub-mode index fastest). See packConformation.

Y = permute(reshape(v, 2, Ns, Nm), [3 2 1]);
Az = Y(:,:,1);
Ar = Y(:,:,2);
end
