function [Gp, Gpp] = linearViscoelasticModuli(pm, freq)
%LINEARVISCOELASTICMODULI Multi-mode Maxwell storage and loss moduli.
%
%   [Gp, Gpp] = linearViscoelasticModuli(pm, freq) returns the polymer
%   contributions of every (MW mode, Zimm sub-mode) pair as
%   (Nm*Ns)-by-numel(freq) matrices; sum over rows for G' and G''.
%   pm comes from modeProperties (fields g, lambda_ms).

[Nm, Ns] = size(pm.lambda_ms);
lamlin = reshape(pm.lambda_ms, [Nm*Ns,1]);
glin = repmat(pm.g, [Ns,1]);
Gp = freq.^2.*glin.*lamlin.^2./(1+(lamlin.*freq).^2);
Gpp = freq.*glin.*lamlin./(1+(lamlin.*freq).^2);
end
