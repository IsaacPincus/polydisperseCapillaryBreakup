# polydisperseCapillaryBreakup

Code for the figures of the paper on capillary thinning of polydisperse
polymer solutions (PS7 and PS16 in DOP). A 0D filament-thinning model with
multi-mode FENE-P polymer stresses is solved for blends whose molecular
weight distributions are taken from GPC data.

## Reproducing a figure

Run the matching script in `PaperScripts/` (from any folder):

| Script                  | Figures      | Run time (desktop) |
|-------------------------|--------------|--------------------|
| `createFig1_FigA1.m`    | 1, A1        | seconds            |
| `createFig2.m`          | 2a, 2b       | seconds            |
| `createFig3.m`          | 3            | < 1 min            |
| `createFig4_Fig5.m`     | 4a, 4b, 5    | minutes            |
| `createFig6a.m`         | 6a           | minutes            |
| `createFig6b.m`         | 6b           | minutes            |
| `createFig7.m`          | 7a-d         | minutes            |
| `createFig8.m`          | 8a, 8b       | ~20 min            |
| `createFig9_Fig10.m`    | 9, 10a, 10b  | minutes            |
| `createFigB1.m`         | B1 (Appendix B) | seconds         |

Each script is three lines: `results = computeFigN(); plotFigN(results);`.
The numbers behind a figure are in `results`. Save the current figure with
`saveFigure('Fig6a')`, which writes .fig/.png/.eps/.emf/.pdf to
`PaperScripts/Figures/`. Some scripts also draw supplementary figures that
are not in the paper; the help text of each `plotFig*` says which is which.

Requires the Curve Fitting Toolbox (`fit`, `prepareCurveData`).

## Layout

```
setupPaths.m                 adds the folders below to the path
Data/                        GPC distributions, SAOS and CaBER measurements
Functions/
  parameters/                defaultParameters (all shared constants), inertialDensity
  data/                      loading the MW distributions and CaBER curves
  polymer/                   MW-distribution splitting, mode moduli and
                             relaxation times, zetacal/ubcal drag model, SAOS moduli
  thinning/                  solveFilamentThinning (the model), relaxation-time
                             extraction for simulations and experiments,
                             concentration sweeps, plate retraction
  util/                      option parsing
  unused/                    functions no figure uses (see below)
PaperScripts/
  createFig*.m               entry points
  compute/computeFig*.m      all calculations for a figure -> results struct
  plot/plotFig*.m            all plotting for a figure, from the results struct
tests/                       regression tests against the original code
```

### Data flow

```
loadMWDistributions --> twoSpeciesModes (one mode per species)  --> modeProperties --> solveFilamentThinning --> simulatedRelaxationTime
                    \-> blendModes (blend split into Nm modes) -/   (g, L, lambda)      (a, Az, Ar, ep, stress)   (t_c, tau_EC)
readThinningCurves --> experimentalRelaxationTime
```

## Parameters

All model parameters are those of Table I of the paper and live in
`Functions/parameters/defaultParameters.m`. Every figure uses them
unchanged; each `computeFig*` only chooses numerical settings (number of
Zimm modes Ns, concentrations, sweep ranges) at its top. Figs 6-8 use the
conformation-dependent drag of Sec. II D; the others constant drag.

The code follows the paper's equations: Eq. 8 (stress balance), Eq. 9
(radius), Eqs. 17-20 (Zimm spectrum with U = 2.36), Eq. 22 (FENE-PM factor
with 3Ns), Eq. 23 (extensibility), Eqs. 26-27 and 31, and Eq. 32 with the
Fig. 3 caption (tau_EC fitted where the strain rate is within 10% of its
minimum).

## Tests

```matlab
cd tests
generateReferenceData()   % once (as slow as running every figure)
runTests()                % function and solver-physics tests (fast)
runTests('Fig3')          % one figure
runTests('all')           % everything
```

* `FunctionTest` compares the functions with the original code (git
  commit `24f47bb`, extracted with `git archive`) on fixed inputs, wherever
  the paper's equations and the original agree (e.g. the solver with Ns = 1).
* `SolverPhysicsTest` checks the solver against the paper's equations
  (f = 1 at equilibrium, Eq. 8 holds at every step, drag ratio >= 1, ...).
* `PaperResultTest` reruns each figure and compares all its numbers with
  those recorded by `generateReferenceData` from the current code, so any
  unintended change to a paper result is caught. After a deliberate change,
  regenerate that figure: `generateReferenceData('Figures', {'Fig3'}, 'Functions', false)`.

Comparison tolerance (`tests/helpers/testTolerance.m`): 1e-9 relative.

## Changes from the original scripts

The refactor itself was checked to reproduce the original scripts bit for
bit. The following were then changed so that the code matches the paper:

* Table I values in every figure: h* = 0.26 (was 0.26, 0.267 or 0.2671),
  a0 = 0.62 mm (was 0.63 mm except Fig 6a), rho = 985 kg/m^3 (was 1000 in
  Figs 3, 6b, 7-9), L from Eq. 23 (Figs 2 and 8 used Lmax/R0).
* Zimm time divided by U = 2.36 (Eq. 17) in Figs 3, 8 and 9 too.
* Eq. 22: f = L^2/(L^2 + 3Ns - sum tr A); the original used L^2 + 3, so
  f was not 1 at equilibrium for Ns > 1.
* Conformation-dependent drag: the ODE used sum(Az - 2Ar) in f; it now uses
  the trace Az + 2Ar like the rest of the model.
* N_K = (L^2 + 1)/3 for the drag model (Sec. II D); was Lmax^2/R0^2.
* tau_EC fit window: 10% of the strain-rate minimum (Fig. 3 caption); was
  20% (Figs 3-6, 9) or 100% (Figs 7c,d, 8).
* Fig 7c,d emphasise the mode whose Wi at the strain-rate minimum just
  exceeds 0.5 (caption); was a fixed mode index (15 or 16).
* Fig 5 has a dashed line at tau_EC = tau_Z and an inset of
  f_{i,min} tau_EC/tau_Z (f_{i,min}: Eq. 22 for species i at the strain-rate
  minimum); Fig 8 marks the sigma_p-weighted curves with symbols, calls the
  Eq. 33 weighting 'terminal', and scales its inset to the simulation.
* The original Fig 9 is now Figs 9 (a0 sweep, Bo > 1 shaded) and 10a (V0
  sweep), with equal y-limits and dotted tau_Z lines; Fig 10b is new: the
  polymer stress of PS16 and PS7 when the plates stop, divided by V0/a0.
  The plate retraction now stops at the gap where the radius equals a0
  (L_end = 8.18 mm for a0 = 0.62 mm), not at a fixed 1.3 mm.
* Plate retraction: the original took the radial velocity of Spiegelberg et
  al. (1996), Eq. 25, which has 1/Lambda0^2 where continuity requires
  1/Lambda0. The strain rate (Eq. B10) and radius (Eq. B9) no longer depend
  on Lambda0 = L0/R0; at the midplane ep = (3/2) (dL/dt)/L and
  R = R0 H^(-3/4).
* Inertia: the term is (C/4) rho a^2 ep^2 with C = 144 as in Eq. 8; the
  original used C rho a^2 ep^2 (4 times larger).
* The strain rate is now found from the stress balance (Eq. 8) in closed
  form: it is a quadratic in ep, and the root that tends to the
  inertialess solution is taken. The original used fzero, which could
  converge to the other, spurious negative root.

Conformation-dependent drag uses Prabhakar's `zetacal` (himodel 2).
Eq. 28 is its large-stretch form: the code uses the full aspect-ratio
correction 1 + ln J(N_t) in place of ln(2 N_t), so that the drag ratio is
1 at equilibrium, caps N_t at N_K and never lets the ratio fall below 1.

## Unused functions

`Functions/unused/` holds functions no paper script calls:
`filamentThinningFENEPM_inertia_initialEXP`, `filamentThinningFENEPM_inertia_simple_CDD`,
`fitGaussiansToMWDist`, `getInitialStretch`, `getStretchExperiments`,
`getStretchExperimentsStressBalanced`, `getTimes`. `tests/UnusedFunctionTest.m`
guards them; delete both together if they are not needed.
