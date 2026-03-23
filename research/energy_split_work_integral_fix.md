# Energy-Split Elastic Energy: Work-Integral Fix

This document records the investigation into a systematic discrepancy
between the analytic energy-split solver and hydro simulation results,
from initial observation through root cause analysis to the implemented
fix.

---

## Table of Contents

1. [Problem Statement](#1-problem-statement)
2. [Initial Observations](#2-initial-observations)
3. [Root Cause Analysis](#3-root-cause-analysis)
4. [Mathematical Derivation](#4-mathematical-derivation)
5. [Two Inconsistencies](#5-two-inconsistencies)
6. [The Fix](#6-the-fix)
7. [Before/After Comparison](#7-beforeafter-comparison)
8. [Conclusions](#8-conclusions)

---

## Test Material Parameters

All tests use a fictitious high-strength material (CGS units):

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Reference density | $\rho_0$ | 4 g/cm$^3$ |
| Bulk sound speed | $C_0$ | 7.5 $\times 10^4$ cm/s |
| Hugoniot slope | $s$ | 1.3 |
| Gruneisen parameter | $\Gamma_0$ | 2.0 |
| Shear modulus | $G$ | 2 $\times 10^{11}$ dyn/cm$^2$ |

Two cases with different yield strength and piston velocity:

| Case | $Y_0$ (dyn/cm$^2$) | $Y_0/(2G)$ | $v_\text{piston}$ (cm/s) | Regime |
|------|---------------------|------------|--------------------------|--------|
| MG Elastic | $1.6 \times 10^{11}$ | 0.40 | 80 000 | Elastic wave |
| MG Two-wave | $1.0 \times 10^{11}$ | 0.25 | 90 000 | Plastic shock + elastic wave |

These have unusually high $Y_0/(2G)$ ratios (0.25--0.40), making the
elastic energy fraction large enough (60--70% of total) to expose the
discrepancy clearly.

---

## 1. Problem Statement

The user reported that their hydro simulation (using energy-split EOS)
does **not converge** to the analytic solver's energy-split solution when
increasing resolution. The discrepancy was visible in both wave positions
and stress/velocity levels:

- **MG Elastic case**: simulation elastic precursor ~2% slower than solver
- **MG Two-wave case**: simulation elastic precursor ~2% slower and
  plastic shock ~4% slower than solver

The user believed the problem was in the analytic solver, not the
simulation.

---

## 2. Initial Observations

Running the solver with `energy_split=True` and comparing against the
simulation images (stress/velocity profiles at $t = 4\;\mu\text{s}$):

**MG Elastic** (stress plot):

| Quantity | Split solver | Simulation (approx.) |
|----------|-------------|---------------------|
| $U_{se}$ | 3095 m/s | ~3025 m/s |
| $\sigma$ | $-0.099$ Mbar | $\sim-0.097$ Mbar |
| Precursor position | 0.01238 m | ~0.0121 m |

**MG Two-wave** (velocity plot):

| Quantity | Split solver | Simulation (approx.) |
|----------|-------------|---------------------|
| $U_{se}$ | 3000 m/s | ~2950 m/s |
| $v_Y$ | 664 m/s | ~650 m/s |
| $U_s$ | 2943 m/s | ~2825 m/s |

The solver systematically overestimates wave speeds. The no-split solver
gives even larger values (3462, 3411 m/s), so the issue is specific to
the energy-split formulation, not a general solver problem.

---

## 3. Root Cause Analysis

The simulation evolves thermal energy directly via a PDE (following
Shashkov 2013):

$$\rho\frac{de_{th}}{dt} = \rho\mathcal{P}_\text{plast} - P\,\nabla\cdot\vec{u}$$

where the plastic dissipation is
$\rho\mathcal{P}_\text{plast} = S_\text{eq}\,\dot\varepsilon^\text{pl}$.
The elastic energy $e_{el}$ is reconstructed diagnostically as
$e_{el} = e - e_{th}$.

The analytic solver instead used the **elastic energy state function**:

$$e_{el}^\text{SF} = \frac{\boldsymbol{S}:\boldsymbol{S}}{4\rho G}$$

and computed $e_{th} = e_\text{total} - e_{el}^\text{SF}$.

These two definitions agree for infinitesimal strains but **diverge for
finite strains**. The discrepancy arises because the time derivative of
the state function and the elastic deviatoric power are not equal:

$$\frac{d}{dt}\left[\frac{\boldsymbol{S}:\boldsymbol{S}}{4\rho G}\right]
= \frac{\boldsymbol{S}:\boldsymbol{\bar D}^{el}}{\rho}
+ \frac{3S_x^2}{8G\rho^2}\frac{d\rho}{dt}$$

The extra term $\frac{3S_x^2}{8G\rho^2}\frac{d\rho}{dt}$ is present in
the state function rate but absent from the elastic power. By using the
state function, the solver attributes this extra energy to the elastic
component, **underestimating** $e_{el}$ and **overestimating** $e_{th}$.
Higher thermal energy means higher EOS pressure, which in turn produces
faster wave speeds.

### Relation to the paper's Section 3.1 identity

The paper's Section 3.1 identifies $e_{el} = \boldsymbol{S}:\boldsymbol{S}/(4\rho G)$
by rewriting the deviatoric work as:

$$\frac{1}{\rho}\,\boldsymbol{S}:\boldsymbol{D}
= \frac{d}{dt}\!\left(\frac{\boldsymbol{S}:\boldsymbol{S}}{4\rho G}\right)
+ \text{plastic term}$$

This step is exact only if $d(\rho G)/dt = 0$, i.e. when the product
of density and shear modulus is constant. Without this approximation,
the exact identity is:

$$\frac{1}{\rho}\,\boldsymbol{S}:\boldsymbol{D}
= \frac{d}{dt}\!\left(\frac{\boldsymbol{S}:\boldsymbol{S}}{4\rho G}\right)
+ \text{plastic term}
+ \frac{\boldsymbol{S}:\boldsymbol{S}}{4\rho G}\,\frac{d}{dt}\ln(\rho G)$$

For the 1D piston problem, $G$ is constant and $\rho$ varies, so the
correction term is:

$$\frac{\boldsymbol{S}:\boldsymbol{S}}{4\rho G}\,\frac{\dot\rho}{\rho}
= \frac{3S_x^2}{8G\rho^2}\,\dot\rho$$

This is precisely the extra term identified above. The paper's
identification of $\boldsymbol{S}:\boldsymbol{S}/(4\rho G)$ as the
elastic energy is therefore a **small-compression approximation**,
whereas the work integral is the exact finite-compression result.

---

## 4. Mathematical Derivation

The correct elastic energy is the integral of the elastic deviatoric
power per unit mass. In the elastic regime, $\boldsymbol{D}^{el} =
\boldsymbol{\bar D}$ (all deformation is elastic), and the elastic power
per unit mass is:

$$\frac{\boldsymbol{S}:\boldsymbol{D}^{el}}{\rho}
= \frac{S_x \cdot \partial v/\partial x}{\rho}
= -\frac{S_x}{\rho^2}\frac{d\rho}{dt}$$

Using the hypoelastic constitutive law $S_x = \frac{4}{3}G\ln(\rho_0/\rho)$,
the elastic energy accumulated from $\rho_0$ to $\rho$ is:

$$e_{el}(\rho) = \int_{\rho_0}^{\rho} \frac{(-S_x)}{\rho'^2}\,d\rho'
= \int_{\rho_0}^{\rho} \frac{4G}{3\rho'^2}\ln\frac{\rho'}{\rho_0}\,d\rho'$$

Substituting $a' = \ln(\rho'/\rho_0)$, $d\rho' = \rho'\,da'$:

$$e_{el} = \frac{4G}{3}\int_0^a \frac{a'}{\rho_0 e^{a'}}\,da'
= \frac{4G}{3\rho_0}\int_0^a a'\,e^{-a'}\,da'$$

Integrating by parts:

$$\boxed{e_{el}(\rho) = \frac{4G}{3\rho_0}\left[1 - (1+a)\,e^{-a}\right],
\qquad a = \ln\frac{\rho}{\rho_0}}$$

### Taylor expansion comparison

At yield, $a_Y = Y_0/(2G)$. The two formulas can be written in
comparable form:

$$e_{el}^{\text{work}} = \frac{4G}{3\rho_0}\left[1 - (1+a_Y)\,e^{-a_Y}\right],
\qquad
e_{el}^{\text{state}} = \frac{2G}{3\rho_0}\,a_Y^2\,e^{-a_Y}$$

Expanding both to cubic order in $a_Y$:

$$e_{el}^{\text{work}} = \frac{2G}{3\rho_0}\,a_Y^2
- \frac{4G}{9\rho_0}\,a_Y^3 + O(a_Y^4)$$

$$e_{el}^{\text{state}} = \frac{2G}{3\rho_0}\,a_Y^2
- \frac{2G}{3\rho_0}\,a_Y^3 + O(a_Y^4)$$

The two agree at $O(a_Y^2)$ but differ at $O(a_Y^3)$: the cubic
coefficients are $-4/9$ (work) vs $-2/3$ (state). The state function
has a larger negative cubic term, making it smaller than the work
integral for finite compression.

Their ratio at leading order is:

$$\frac{e_{el}^\text{work}}{e_{el}^\text{state}} \approx \frac{\rho}{\rho_0}
= e^{a_Y} = e^{Y_0/(2G)}$$

For $Y_0/(2G) = 0.25$ (the two-wave case): $e^{0.25} = 1.284$, so the
work integral exceeds the state function by ~28% at leading order.
The actual ratio at yield is:

- State function: $e_{el}^Y = (2G/3\rho_0)\,a_Y^2\,e^{-a_Y} = 1.622 \times 10^9$
- Work integral: $e_{el}^Y = 1.767 \times 10^9$ (9% larger)

The leading-order ratio overestimates the correction because the cubic
and higher terms partly cancel the $e^{a_Y}$ factor.

---

## 5. Two Inconsistencies

### 5a. Elastic regime

In the elastic regime ($P_\text{plast} = 0$), the simulation's thermal
energy PDE reduces to:

$$\rho\frac{de_{th}}{dt} = -P\,\nabla\cdot\vec{u} = \frac{P}{\rho}\frac{d\rho}{dt}$$

Only PdV work heats the material. The elastic energy is the work
integral above. But the solver's state function gives $e_{th} =
e_\text{total} - S^2/(4\rho G)$, which assigns **less** to elastic and
**more** to thermal, inflating the EOS pressure and wave speed.

### 5b. Plastic regime

At yield, the deviatoric stress is constant ($S_x = -\frac{2}{3}Y_0$)
and all deformation is plastic ($\boldsymbol{D}^{el} = 0$). Therefore:

$$\frac{de_{el}}{dt} = \frac{\boldsymbol{S}:\boldsymbol{D}^{el}}{\rho} = 0$$

The elastic energy is **frozen** at its yield value. The thermal energy
PDE gives $de_{th} = de_\text{total}$ (plastic dissipation equals total
deviatoric power), so $e_{el} = e - e_{th}$ does not change across the
plastic shock.

But the old solver recomputed:

$$e_{el,2}^\text{old} = \frac{Y_0^2}{6\rho_2 G}$$

Since $\rho_2 > \rho_Y$, this gives $e_{el,2} < e_{el}^Y$, incorrectly
reducing the elastic energy and inflating $e_{th,2}$.

### Physical interpretation

Freezing $e_{el}$ does not violate energy conservation. The quantity
$e_{el}$ is a *specific* energy (per unit mass). For a material parcel
of fixed mass $m$, the total elastic energy is:

$$E_{el}^\text{parcel} = m\,e_{el}$$

If $e_{el}$ is frozen after yield, then the total deviatoric elastic
energy of that parcel also remains constant:

$$\frac{dE_{el}^\text{parcel}}{dt} = m\,\frac{de_{el}}{dt} = 0$$

Further compression still increases the total internal energy of the
parcel, but in this scheme that increase is attributed entirely to
pressure work and plastic dissipation, not to additional reversible
deviatoric storage. The elastic "bucket" is full at yield and does not
grow further.

---

## 6. The Fix

### 6a. New helper method

Added `_elastic_energy_work_integral(self, rho)` computing:

$$e_{el}(\rho) = \frac{4G}{3\rho_0}\left[1 - (1+a)\,e^{-a}\right]$$

### 6b. Elastic regime changes

- `_solve_yield_and_precursor`: replaced
  `e_el_Y = Y_0**2 / (6*rho_Y*G)` with
  `e_el_Y = self._elastic_energy_work_integral(rho_Y)`
- `_solve_elastic_wave` residual: replaced
  `e_eos = e_E - 1.5*Sx**2/(4*rho_E*G)` with
  `e_eos = e_E - self._elastic_energy_work_integral(rho_E)`
- Post-solve `e_el_Y` assignment: same work-integral call

### 6c. Plastic regime changes

- `f_shock_02` and `f_shock_y2` residuals: replaced
  `e_2 - Y_0**2/(6*rho_2*G)` with `e_2 - self.e_el_Y` (frozen at yield)
- `_solve_shocked_state_y2` and `_solve_shocked_state_02`: replaced
  `e_el_2 = Y_0**2/(6*rho_2*G)` with `e_el_2 = self.e_el_Y`

### 6d. No changes to non-split mode

All modifications are gated behind `if self.energy_split`. The
`energy_split=False` path is completely unchanged.

---

## 7. Before/After Comparison

Simulation reference values are approximate, extracted from user-provided
plots at $t = 4\;\mu\text{s}$.

### MG Elastic ($Y_0 = 1.6 \times 10^{11}$, $v_p = 80\,000$ cm/s)

| Quantity | Old solver | New solver | Simulation | Old err | New err |
|----------|-----------|-----------|------------|---------|---------|
| $U_{se}$ (cm/s) | 309 454 | 304 177 | ~302 500 | +2.3% | +0.6% |
| Precursor at 4 $\mu$s (m) | 0.01238 | 0.01217 | ~0.0121 | +2.3% | +0.6% |
| $\sigma$ (Mbar) | $-0.0990$ | $-0.0973$ | $\sim-0.097$ | +2.1% | +0.3% |
| $e_{el}^Y$ (erg/g) | $2.211 \times 10^9$ | $2.540 \times 10^9$ | — | — | — |
| $e_{th}^Y$ (erg/g) | $9.888 \times 10^8$ | $6.605 \times 10^8$ | — | — | — |

### MG Two-wave ($Y_0 = 1.0 \times 10^{11}$, $v_p = 90\,000$ cm/s)

| Quantity | Old solver | New solver | Simulation | Old err | New err |
|----------|-----------|-----------|------------|---------|---------|
| $U_{se}$ (cm/s) | 300 016 | 296 097 | ~295 000 | +1.7% | +0.4% |
| $v_Y$ (cm/s) | 66 363 | 65 496 | ~65 000 | +2.1% | +0.8% |
| $U_s$ (cm/s) | 294 308 | 283 028 | ~282 500 | +4.2% | +0.2% |
| Precursor at 4 $\mu$s (m) | 0.01200 | 0.01184 | ~0.0118 | +1.7% | +0.4% |
| Shock at 4 $\mu$s (m) | 0.01177 | 0.01132 | ~0.0113 | +4.2% | +0.2% |
| $P_Y$ (dyn/cm$^2$) | $1.297 \times 10^{10}$ | $1.091 \times 10^{10}$ | — | — | — |
| $e_{el}^Y$ (erg/g) | $1.623 \times 10^9$ | $1.767 \times 10^9$ | — | — | — |
| $e_{el,2}$ (erg/g) | $1.454 \times 10^9$ | $1.767 \times 10^9$ | — | — | — |

### Key improvements

| Quantity | Old error | New error |
|----------|----------|----------|
| MG Elastic $U_{se}$ | +2.3% | **+0.6%** |
| MG Elastic $\sigma$ | +2.1% | **+0.3%** |
| MG Two-wave $U_{se}$ | +1.7% | **+0.4%** |
| MG Two-wave $v_Y$ | +2.1% | **+0.8%** |
| MG Two-wave $U_s$ | +4.2% | **+0.2%** |

The largest improvement is in the plastic shock speed ($U_s$), which
went from 4.2% error to 0.2%. The remaining ~0.3--0.8% residual is
within the uncertainty of reading wave positions from the simulation
plots.

### Impact on typical metals

For typical metals ($Y_0/(2G) \sim 10^{-3}$), the correction factor
$e^{Y_0/(2G)} - 1 \approx Y_0/(2G) \sim 0.001$, so the state function
and work integral agree to $< 0.1\%$. The fix is only significant for
materials with $Y_0/(2G) \gtrsim 0.05$.

---

## 8. Conclusions

1. **The elastic energy state function $\boldsymbol{S}:\boldsymbol{S}/(4\rho G)$
   is inconsistent with the thermal energy PDE used in hydrocodes.** The
   correct elastic energy is the work integral of the elastic deviatoric
   power, which differs from the state function by a density-dependent
   correction that grows with the elastic strain $a = \ln(\rho/\rho_0)$.

2. **During plastic deformation, the elastic energy is frozen at its
   yield-point value**, because $\boldsymbol{D}^{el} = 0$ when the
   deviatoric stress is constant at yield. The old solver incorrectly
   recomputed $e_{el} = Y_0^2/(6\rho_2 G)$ with the post-shock density.

3. **The fix resolves the simulation-solver discrepancy.** Errors dropped
   from 1.7--4.2% to 0.2--0.8%, consistent with the image-reading
   uncertainty.

4. **The correction is negligible for typical metals** ($Y_0/(2G) \ll 1$)
   and only matters for high-strength or high-strain materials where
   $Y_0/(2G) \gtrsim 0.05$.

5. **The non-split mode is unaffected.** All changes are gated behind
   `energy_split=True`.

---

## Relation to Previous Investigations

This fix addresses a different issue from the earlier root-finding
problem documented in
[`research/Opus/summary.md`](Opus/summary.md) and
[`research/Codex/summary.md`](Codex/summary.md). Those investigations
focused on the **root landscape** of $f_\text{shock}(U_s)$ -- cases
where the weak-shock root vanished and the solver grabbed an unphysical
strong-shock root. The fix there was robust multi-root scanning.

The current issue is about the **energy decomposition formula** itself:
even when the correct root is found, the solver produced slightly wrong
values because $e_{el}$ was computed from the state function rather than
the work integral. The two fixes are complementary and independent.

---

## Files Modified

| File | Changes |
|------|---------|
| `elastoplastic_piston_solver.py` | Added `_elastic_energy_work_integral`; replaced state function with work integral in yield solve, elastic wave solve, plastic shock residuals, and shocked state assignments |
| `README.md` | Updated energy-split formulas and documentation |
