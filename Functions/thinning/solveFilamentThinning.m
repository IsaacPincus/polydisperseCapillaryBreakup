function sol = solveFilamentThinning(p, modes, varargin)
%SOLVEFILAMENTTHINNING 0D capillary thinning of a multi-mode FENE-P filament.
%
%   sol = solveFilamentThinning(p, modes, Name, Value)
%
%   Integrates the filament radius a(t) and the conformation tensor
%   components Az, Ar of every mode, with the strain rate ep found at each
%   instant from the stress balance
%       rho a^2 ep^2/8 + sum(g f (Az - Ar)) + 3 eta_s ep = (2X-1) chi/a,
%   until a/a0 falls to p.cutoff or t reaches p.tmax.
%
%   Inputs
%     p      parameters: chi, X, eta_s, rho, a0, cutoff, tmax (and hsk,
%            phieq for conformation-dependent drag)
%     modes  from modeProperties: g (Nm-by-1), L (Nm-by-1), lambda_ms
%            (Nm-by-Ns), and Nk (Nm-by-1) for conformation-dependent drag
%
%   Options
%     'Drag'               'constant' (default) or 'conformationDependent'.
%                          The latter scales the slowest sub-mode of each
%                          species by zetacal(..., himodel = 2), Sec. II D.
%     'InitialConformation' state vector of length 2*Nm*Ns ordered
%                          [Az(1,1) Ar(1,1) Az(1,2) Ar(1,2) ... Az(Nm,Ns) Ar(Nm,Ns)]
%                          (default: all ones, i.e. equilibrium)
%     'InitialStrainRate'  initial guess for ep at t = 0
%                          (default: sqrt(8 chi/(a0^3 rho)))
%
%   Output fields
%     t, a        time and radius (column vectors)
%     Az, Ar      length(t)-by-Nm-by-Ns
%     f           FENE-P factor, length(t)-by-Nm
%     stress      polymer stress of each MW mode, length(t)-by-Nm
%     ep          strain rate (column)
%     Wi          effective Weissenberg number of the slowest sub-mode, Nm-by-length(t)
%     te, ye, ie  event output of ode15s
%
%   Equation numbers refer to the paper: the stress balance is Eq. 8, the
%   radius Eq. 9, the conformations Eqs. 26 and 31 with the FENE-PM factor
%   f = L^2/(L^2 + 3 Ns - sum_s tr A_s) of Eq. 22, and the stress Eq. 27.

[Nm, Ns] = size(modes.lambda_ms);
opt = parseOptions(struct('Drag', 'constant', ...
    'InitialConformation', ones(2*Nm*Ns, 1), ...
    'InitialStrainRate', []), varargin{:});

c.Nm = Nm;
c.Ns = Ns;
c.chi = p.chi;
c.mu = p.eta_s;
c.X = p.X;
c.rho = p.rho;
c.a0 = p.a0;
c.L2 = modes.L.^2;
c.trace0 = 3*Ns;   % sum of tr(A_s) at equilibrium, Eq. 22
c.g = modes.g;
c.lambda = modes.lambda_ms;
switch opt.Drag
    case 'constant'
        c.cdd = false;
    case 'conformationDependent'
        c.cdd = true;
        c.Nk = modes.Nk;
        c.hsk = p.hsk;
        c.phieq = p.phieq;
    otherwise
        error('solveFilamentThinning:drag', 'Unknown Drag "%s".', opt.Drag);
end
if isempty(opt.InitialStrainRate)
    c.epReset = sqrt(8*c.chi/(c.a0^3*c.rho));
else
    c.epReset = opt.InitialStrainRate;
end

y0 = [c.a0; opt.InitialConformation(:)];
options = odeset('Events', @(t,y) radiusDecayed(t, y, c.a0, p.cutoff), ...
    'RelTol', 1e-6, 'AbsTol', 1e-9);
[t, y, te, ye, ie] = ode15s(@(t,y) thinningRHS(t, y, c), [0, p.tmax], y0, options);

sol = postProcess(t, y, c);
sol.te = te;
sol.ye = ye;
sol.ie = ie;
end

% -------------------------------------------------------------------------
function dydt = thinningRHS(~, y, c)
persistent epGuess
a = y(1);
[Az, Ar] = unpackConformation(y(2:end), c.Nm, c.Ns);

f = feneFactor(Az, Ar, c);
if c.cdd
    zet = dragRatio(Az, Ar, c);
else
    zet = ones(c.Nm, c.Ns);
end

% fzero is seeded with the previous strain rate; reset at the start of a run
if a == c.a0
    epGuess = c.epReset;
end
P = sum(c.g.*f.*sum(Az-Ar,2));   % Eq. 27
ep = fzero(@(ep) stressBalance(ep, a, P, c), epGuess);
epGuess = ep;

adot = -0.5*ep*a;
Azdot = 2*ep*Az - f./(c.lambda.*zet).*(Az-1);
Ardot = -ep*Ar - f./(c.lambda.*zet).*(Ar-1);
dydt = [adot; packConformation(Azdot, Ardot)];
end

function sol = postProcess(t, y, c)
lent = length(t);
Y = reshape(y(:,2:end), lent, 2, c.Ns, c.Nm);
Az = permute(Y(:,1,:,:), [1 4 3 2]);
Ar = permute(Y(:,2,:,:), [1 4 3 2]);
a = y(:,1);

L2 = c.L2;
f = L2'./(L2'+c.trace0-sum(Az+2*Ar,3));
stress = c.g'.*f.*sum(Az-Ar,3);

ep = zeros(size(t));
guess = sqrt(8*c.chi/(a(1)^3*c.rho));
for tt = 1:lent
    Azt = reshape(Az(tt,:,:), c.Nm, c.Ns);
    Art = reshape(Ar(tt,:,:), c.Nm, c.Ns);
    P = sum(c.g.*feneFactor(Azt, Art, c).*sum(Azt-Art,2));
    ep(tt) = fzero(@(ep) stressBalance(ep, a(tt), P, c), guess);
    guess = ep(tt);
end

Wi = c.lambda(:,1).*ep';
if c.cdd
    zet = zeros(lent, c.Nm);
    for mm = 1:c.Nm
        Q = sqrt((Az(:,mm,1) + 2*Ar(:,mm,1))/3);
        for tt = 1:lent
            zet(tt,mm) = zetacal(Q(tt), 0, c.phieq, c.Nk(mm), c.hsk, 2);
        end
    end
    Wi = Wi.*zet';
end

sol = struct('t', t, 'a', a, 'Az', Az, 'Ar', Ar, 'f', f, ...
    'stress', stress, 'ep', ep, 'Wi', Wi);
end

function f = feneFactor(Az, Ar, c)
% FENE-PM factor of each species, Eq. 22
f = c.L2./(c.L2+c.trace0-sum(Az+2*Ar,2));
end

function bal = stressBalance(ep, a, P, c)
% Eq. 8 with c.rho = C*rho*8 (see inertialDensity)
bal = 1/8*ep^2*a^2*c.rho + P ...
    + 3*c.mu*ep - (2*c.X-1)*c.chi/a;
end

function zet = dragRatio(Az, Ar, c)
% conformation-dependent drag of the slowest sub-mode of each species
zet = ones(c.Nm, c.Ns);
Q = sqrt((Az(:,1) + 2*Ar(:,1))/3);
for mm = 1:c.Nm
    zet(mm,1) = zetacal(Q(mm), 0, c.phieq, c.Nk(mm), c.hsk, 2);
end
end

function [value, isterminal, direction] = radiusDecayed(~, y, a0, cutoff)
value = y(1)/a0 - cutoff;
isterminal = 1;
direction = 0;
end
