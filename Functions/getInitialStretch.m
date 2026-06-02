function [Az,Ar] = getInitialStretch(L0,L1,R0)
    %GETINITIALSTRETCH Calculates the affine stretching from the initial
    %plate retraction in a CaBER experiment, based upon a cylindrical
    %intiial shape for the fluid and a catenoidal final shape. This may
    %fail at Bo>0!

    a0 = R0;
    a = fzero(@(a) catenoidEquation(a,L0,L1,R0), a0);
    R1 = R0 - a*cosh(L1/(2*a));
    
    % assume Az and Ar are initially 1
    Az = (L1/L0)^2;
    Ar = (R1/R0)^2;
end

function f = catenoidEquation(a, L0,L1,R0)
    R1 = R0 - a*cosh(L1/(2*a));
    f = pi*R0^2*L0 ...
        - 1/2*(a^3*sinh(L1/a) + L1*(a^2+2*R1^2) + 8*a^2*R1*sinh(L1/(2*a)));
end

