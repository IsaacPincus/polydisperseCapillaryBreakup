function [zet, G, zb, Nb] = zetacal(Q, ~, phieq, nk, hsk, himodel)
%ZETACAL Conformation-dependent friction coefficient of a chain.
%
%   Written by R. Prabhakar (Prabhakar Ranganathan), who provided the
%   conformation-dependent drag implementation used in the paper (Sec. II D).
%   See R. Prabhakar et al., J. Rheol. 60, 345-366 (2016), and
%   R. Prabhakar et al., Phys. Rev. Fluids 2, 011301 (2017).
%   Restructured (shared blocks moved into local functions) without
%   changing any results.
%
%   [zet, G, zb, Nb] = zetacal(Q, d, phieq, nk, hsk, himodel)
%
%   Outputs
%     zet  friction coefficient scaled by the Zimm friction coefficient
%          (isolated, isotropic coil at equilibrium)
%     G    correction to the Rouse-like drag of a chain of blobs
%     zb   friction coefficient of the inner-most blob
%     Nb   number of inner-most blobs
%
%   Inputs
%     Q       rms end-to-end distance scaled by its equilibrium value Q0
%     d       transverse coil diameter (unused; kept for the call signature)
%     phieq   equilibrium coil volume fraction c/c*
%     nk      number of Kuhn segments
%     hsk     hydrodynamic interaction parameter of a Kuhn segment
%     himodel 0 - constant (freely draining Rouse) drag, no HI
%             1 - constant drag of a fully stretched chain, with
%                 concentration-dependent Batchelor correction
%             2 - conformation-dependent drag, no concentration dependence
%                 (the model used for the paper)
%             3 - conformation- and concentration-dependent drag

if hsk == 0
    hsk = eps;
end
L = nk^0.5;

switch himodel
    case 0
        zet = 2*pi^0.5*hsk*nk^0.5/ubcal(hsk, nk);
        G = 1;  zb = 1;  Nb = 1;   % not meaningful for this model

    case 1
        [zeta_Rouse, G, zb] = fullyStretchedDrag(nk, hsk, phieq);
        zet = zeta_Rouse/G;
        Nb = nk;

    case 2
        zeta_isolcoil = ubcal(hsk, nk);      % equilibrium coil drag / (3 pi eta_s Q0)
        xip = max(1/Q, 1/nk^0.5);            % tension-blob size
        Np = min(Q^2, nk);                   % number of tension blobs
        A = Q/xip;                           % aspect ratio of the blob pole
        zb = ubcal(hsk, nk/Np)*xip;          % blob drag / (3 pi eta_s Q0)
        zb = zb/zeta_isolcoil;               % ... / Zimm friction coefficient
        Nb = Np;
        zeta_Rouse = zb*Np;
        G = interBlobCorrection(A);
        zet = zeta_Rouse/G;
        if zet < 1
            zet = 1;
        end

    case 3
        zeta_isolcoil = ubcal(hsk, nk);
        % coiled-state drag
        if phieq > L        % concentrated (unentangled)
            zb = ubcal(hsk, 1)*1/sqrt(nk)/zeta_isolcoil;
            zeta_coil = zb*nk;
            Nc = nk;
        elseif phieq > 1    % semidilute: chain of concentration blobs
            Nc = phieq^2;
            xic = 1/phieq;
            zb = ubcal(hsk, nk/Nc)*xic;
            zb = zb/zeta_isolcoil;
            zeta_coil = zb*Nc;
        else
            zeta_coil = 1;
            Nc = 1;
            zb = 1;
        end

        if Q < nk^0.5      % partially stretched chain
            xip = 1/Q;
            Np = Q^2;
            A = Q/xip;
            if phieq < Q    % weak screening
                Nb = Np;
                zb = ubcal(hsk, nk/Np)*xip;
                zb = zb/zeta_isolcoil;
                zeta_Rouse = zb*Np;
                if phieq > 0
                    G = interBlobCorrection(A, phieq*(Q*xip^2));
                else
                    G = interBlobCorrection(A);
                end
            else            % strong screening: no inter-blob HI
                Nb = Nc;
                zeta_Rouse = zeta_coil;
                G = 1;
            end
        else               % Gaussian chain stretched past its contour length
            [zeta_Rouse, G, zb] = fullyStretchedDrag(nk, hsk, phieq);
            Nb = nk;
        end
        zet = zeta_Rouse/G;
end

zet = abs(zet);
end

function [zeta_Rouse, G, zb] = fullyStretchedDrag(nk, hsk, phieq)
% Drag of a fully stretched chain (Q = nk^0.5) of nk tension blobs.
xip = 1/nk^0.5;
Np = nk;
A = nk;
Q = nk^0.5;
zb = ubcal(hsk, nk/Np)*xip/ubcal(hsk, nk);
zeta_Rouse = zb*Np;
if phieq > 0 && phieq < nk^0.5
    G = interBlobCorrection(A, phieq*(Q*xip^2));
elseif phieq == 0
    G = interBlobCorrection(A);
else
    G = 1;
end
end

function G = interBlobCorrection(A, phipincus)
% Correction for hydrodynamic interactions between blobs in a rod of
% aspect ratio A. With phipincus (volume fraction of blob rods), screening
% follows Mackaplow & Shaqfeh, JFM 329 (1996) via a Batchelor-style
% interpolation; otherwise G = 1 + log(J).
e1 = exp(1);
J = 1 + (0.5*(A-1) + (2/e1)*(A-1)*(A-1))/A;
if nargin < 2
    G1 = 1 + log(J);
    G2 = 1;
else
    chibyxip = (log(1/phipincus)/phipincus)^0.5;  % screening length / blob size
    B = chibyxip + 1;
    G1 = 1 + log((J*(B-1) + A)/(B + A));
    G2 = 1 + log(((B-1) + 1)/(B + 1/A));          % gives G = 1 in limiting cases
end
G = G1/G2;
end
