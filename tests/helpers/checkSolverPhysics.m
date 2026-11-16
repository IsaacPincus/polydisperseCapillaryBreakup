function problems = checkSolverPhysics()
%CHECKSOLVERPHYSICS Check solveFilamentThinning against the paper's equations.
%
%   problems = checkSolverPhysics() runs small multi-mode problems (two
%   species, Ns = 3) with constant and conformation-dependent drag and
%   returns a message for every violated property (empty if all hold):
%     * Eq. 22: f = 1 at equilibrium for every species
%     * Eq. 8:  the stress balance holds at every output time
%     * Eq. 9:  the radius decreases and the run stops at a/a0 = cutoff
%     * Eq. 31: the drag ratio is >= 1, so Wi with CDD >= tau_Z*ep

problems = {};
p = defaultParameters();
mw = loadMWDistributions(p);
pm = modeProperties(twoSpeciesModes(mw, 1000, 100, p), p, 3);

% Eq. 22 at equilibrium, through the solver's own output at t = 0
for drag = {'constant', 'conformationDependent'}
    sol = solveFilamentThinning(p, pm, 'Drag', drag{1});
    name = drag{1};
    if any(abs(sol.f(1,:) - 1) > 1e-12)
        problems{end+1} = sprintf('%s: f(t=0) = %s, expected 1 (Eq. 22)', name, mat2str(sol.f(1,:))); %#ok<AGROW>
    end

    % Eq. 8 residual, scaled by the capillary pressure
    capillary = (2*p.X-1)*p.chi./sol.a;
    polymer = sum(sol.stress, 2);
    residual = 1/8*sol.ep.^2.*sol.a.^2*p.rho + polymer + 3*p.eta_s*sol.ep - capillary;
    if max(abs(residual)./capillary) > 1e-8
        problems{end+1} = sprintf('%s: stress balance violated by %.2g (Eq. 8)', ...
            name, max(abs(residual)./capillary)); %#ok<AGROW>
    end

    % Eq. 9 and the stopping criterion
    if any(diff(sol.a) > 0)
        problems{end+1} = sprintf('%s: radius increases', name); %#ok<AGROW>
    end
    if isempty(sol.te) || abs(sol.a(end)/(p.a0*p.cutoff) - 1) > 1e-3
        problems{end+1} = sprintf('%s: did not stop at a/a0 = cutoff', name); %#ok<AGROW>
    end

    % stress definition, Eq. 27
    stress = pm.g'.*sol.f.*sum(sol.Az - sol.Ar, 3);
    if max(abs(stress(:) - sol.stress(:))) > 1e-12*max(abs(sol.stress(:)))
        problems{end+1} = sprintf('%s: stress differs from Eq. 27', name); %#ok<AGROW>
    end

    if strcmp(name, 'conformationDependent')
        WiConstant = pm.lambda_ms(:,1).*sol.ep';
        if any(sol.Wi(:) < WiConstant(:)*(1 - 1e-12))
            problems{end+1} = 'conformationDependent: drag ratio below 1 (Eq. 31)'; %#ok<AGROW>
        end
    end
end
end
