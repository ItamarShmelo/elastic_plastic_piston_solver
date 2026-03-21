# Investigation Summary: Unphysical Solutions in the Energy-Split Elastoplastic Piston Solver

This document records the full investigation into an unphysical solution produced
by the energy-split mode of the 1-D elastoplastic piston solver, from initial
diagnosis through solver fix, parameter-space mapping, EOS comparison, and wave
speed analysis.

All scripts that produced the results below are in [`scripts/`](scripts/).
A local copy of the modified solver is included for reproducibility.

---

## Table of Contents

1. [Problem Statement](#1-problem-statement)
2. [Phase 1 -- Root Landscape Diagnosis](#2-phase-1----root-landscape-diagnosis)
3. [Phase 2 -- Solver Fix](#3-phase-2----solver-fix)
4. [Phase 3 -- Parameter Space Exploration](#4-phase-3----parameter-space-exploration)
5. [Phase 4 -- Regime Maps (Mie-Gruneisen EOS)](#5-phase-4----regime-maps-mie-gruneisen-eos)
6. [Phase 5 -- Physical Interpretation of Regimes](#6-phase-5----physical-interpretation-of-regimes)
7. [Phase 6 -- Ideal Gas EOS Regime Maps](#7-phase-6----ideal-gas-eos-regime-maps)
8. [Phase 7 -- Wave Speed Ratio Analysis](#8-phase-7----wave-speed-ratio-analysis)
9. [Phase 8 -- Ideal Gas Profile Comparison](#9-phase-8----ideal-gas-profile-comparison)
10. [Phase 9 -- MG-to-Ideal-Gas Parameter Mapping](#10-phase-9----mg-to-ideal-gas-parameter-mapping)
11. [Conclusions](#11-conclusions)

---

## Material Parameters (CGS)

Unless otherwise noted, all runs use these aluminum-like parameters:

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Reference density | $\rho_0$ | 2.79 g/cm$^3$ |
| Bulk sound speed | $C_0$ | 5.330 $\times 10^5$ cm/s |
| Hugoniot slope | $s$ | 1.34 |
| Gruneisen parameter | $\Gamma_0$ | 2.0 |
| Specific heat | $c_v$ | 9.0 $\times 10^6$ erg/(g K) |
| Shear modulus | $G$ | 2.86 $\times 10^{11}$ dyn/cm$^2$ |
| Yield stress | $Y_0$ | variable ($Y_0/G$ swept) |

---

## 1. Problem Statement

The elastoplastic piston solver supports two modes:

- **Original mode**: the full specific internal energy $e$ is fed to the
  equation of state.
- **Energy-split mode**: the internal energy is decomposed as
  $e = e_\mathrm{th} + e_\mathrm{el}$, where
  $e_\mathrm{el} = S_{ij}S_{ij}/(4\rho G) = Y_0^2/(6\rho G)$,
  and only $e_\mathrm{th}$ is passed to the EOS.

At $Y_0/G = 0.15$ and $v_\mathrm{piston} = 80\,000$ cm/s, the energy-split
mode produced an **unphysical** shocked state: extreme compression
($\rho_2/\rho_Y \approx 2$), a shock speed far below the elastic precursor
($U_s/U_{se} \approx 0.2$), and proximity to the Mie-Gruneisen Hugoniot
singularity ($\mu_2 \to 1/s$).

The hypothesis was that $f_\mathrm{shock}(U_s)$ has **multiple roots**, and the
solver was selecting the wrong one.

---

## 2. Phase 1 -- Root Landscape Diagnosis

**Script**: [`scripts/diagnose_energy_split.py`](scripts/diagnose_energy_split.py)

The diagnostic script evaluates $f_\mathrm{shock}(U_s)$ across the full range
$(v_\mathrm{piston}, U_{se})$ with 10 000 probe points for both the original and
energy-split modes.

### Original mode

![f_shock landscape -- original mode](assets/diag_fshock_orig.png)

In the **original mode**, $f_\mathrm{shock}$ starts positive near $U_{se}$,
crosses zero once (the physical weak-shock root at $U_s/U_{se} \approx 0.97$),
and remains negative toward $v_\mathrm{piston}$. The solver correctly finds the
single root with moderate compression ($\rho_2/\rho_Y \approx 1.04$).

### Energy-split mode

![f_shock landscape -- energy-split mode](assets/diag_fshock_split.png)

In the **energy-split mode**, $f_\mathrm{shock}$ starts **negative** near
$U_{se}$ (the weak-shock root has vanished), then crosses zero once at
$U_s/U_{se} \approx 0.19$, producing only a strong-shock root with:

- $\rho_2/\rho_Y \approx 1.97$ (near 2x compression)
- $\mu_2 \approx 0.65$ (close to the Hugoniot singularity at $1/s = 0.746$)

The old solver's downward search from $U_{se}$ found this strong-shock root and
returned it as the "solution", producing the unphysical state.

### Root analysis at the diagnostic parameters

| Quantity | Original mode | Energy-split mode |
|----------|--------------|-------------------|
| Number of roots | 1 | 1 |
| $U_s$ (cm/s) | $\approx 5.4 \times 10^5$ | $\approx 1.1 \times 10^5$ |
| $U_s/U_{se}$ | $\approx 0.97$ | $\approx 0.19$ |
| $\rho_2/\rho_Y$ | $\approx 1.04$ | $\approx 1.97$ |
| Physical? | Yes (weak shock) | No (strong shock) |

The key insight: the energy-split mode shifts $f_\mathrm{shock}(U_{se}^-)$ from
positive to negative, destroying the weak-shock root. The function still crosses
zero once (the strong-shock branch), and the solver grabs it.

---

## 3. Phase 2 -- Solver Fix

**Modified file**: [`scripts/elastoplastic_piston_solver.py`](scripts/elastoplastic_piston_solver.py)

The root-finding in `_solve_wave_structure` was replaced with a new method
`_find_physical_shock_speed` that:

1. **Scans the full range** $(v_\mathrm{piston}, U_{se})$ with 2000 linear
   probes, locating every sign change of $f_\mathrm{shock}$.
2. **Solves all brackets** with Brent's method.
3. **Filters by physical constraints**:
   - $\rho_2 > \rho_Y$ (monotonic compression)
   - $\mu_2 < 1/s$ (below Hugoniot singularity)
4. **Selects the largest $U_s$** (weakest shock, least compression) -- this is
   the physical root for the piston problem.
5. **Warns** if only the strong-shock root survives (heuristic:
   $U_s < 0.5\,U_{se}$ with only one valid root), indicating an unphysical
   solution.

```python
def _find_physical_shock_speed(self, f_shock, v_piston, v_Y, U_se, rho_Y, n_scan=2000):
    # Dense linear scan to locate every sign change
    probes = np.linspace(v_piston * (1 + 1e-8), U_se * (1 - 1e-10), n_scan)
    # ... find all brackets, solve each, filter by physics, pick max(U_s)
```

This replaced the old downward geometric search that found only the first sign
change from $U_{se}$ -- which worked when $f_\mathrm{shock}(U_{se}^-) > 0$ but
failed when the weak-shock root vanished.

---

## 4. Phase 3 -- Parameter Space Exploration

With the fixed solver, we explored how $Y_0/G$ and $v_\mathrm{piston}$ affect
which root exists:

- **Reducing $Y_0/G$** (at fixed $v_\mathrm{piston} = 80\,000$ cm/s):
  below $Y_0/G \approx 0.12$ the physical weak-shock root reappears in
  energy-split mode.
- **Reducing $v_\mathrm{piston}$** (at fixed $Y_0/G = 0.15$): below
  $v_p \approx 78\,000$ cm/s the weak-shock root reappears.
- The physical 2-wave solution exists in a **triangular wedge** in
  $(v_\mathrm{piston},\, Y_0/G)$ space, bounded below by $v_Y(Y_0/G)$ and
  above by a critical velocity $v_p^*(Y_0/G)$.

---

## 5. Phase 4 -- Regime Maps (Mie-Gruneisen EOS)

**Script**: [`scripts/regime_map.py`](scripts/regime_map.py)

We classified every point on a 2-D grid of $Y_0/G \in [0.005, 0.80]$ and
$v_\mathrm{piston} \in [2000, 400\,000]$ cm/s into one of five regimes:

| # | Regime | Color | Description |
|---|--------|-------|-------------|
| 1 | **No plastic shock** | Blue | $v_p < v_Y$: purely elastic response |
| 2 | **Physical 2-wave** | Teal | Two roots; weak shock selected. Elastic precursor + plastic shock coexist |
| 3 | **Strong shock only** | Yellow | Only the strong-shock root survives (unphysical) |
| 4 | **No root** | Red | $f_\mathrm{shock}$ has no sign change anywhere |
| 5 | **Overdriven** | Purple | $v_p \geq U_{se}$: piston outruns the elastic precursor |

### Energy-split mode

![MG energy-split regime map](assets/regime_map_split.png)

### Original mode

![MG original mode regime map](assets/regime_map_orig.png)

### Combined overlay (energy-split grid, both boundaries)

![MG combined regime map](assets/regime_map_combined.png)

### Regime statistics (Mie-Gruneisen)

| Regime | Energy-split (%) | Original (%) |
|--------|-----------------|--------------|
| No plastic shock | 34.4 | 35.1 |
| Physical 2-wave | 11.7 | 12.5 |
| Strong shock only | 46.2 | 45.2 |
| No root | 7.7 | 7.1 |
| Overdriven | 0.0 | 0.0 |

### Key observations

1. The **physical 2-wave wedge** is narrow: only ~12% of the parameter space.
   It is bounded below by $v_Y(Y_0/G)$ (diagonal) and above by
   $v_p^* \approx 0.163\,C_0 \approx 86\,900$ cm/s (near-horizontal).

2. The **strong-shock-only** (yellow) region dominates, occupying ~46%.
   This is the region where the weak-shock root has vanished due to the
   Mie-Gruneisen EOS behavior.

3. The energy-split mode shifts the 2-wave/strong-shock boundary to
   **lower** $v_\mathrm{piston}$ (by up to ~10% at high $Y_0/G$).

4. The **no-root** (red) region appears at extreme velocities
   ($v_p/C_0 \gtrsim 0.65$) due to the MG Hugoniot singularity.

### Boundary table (excerpt): maximum $v_p$ for physical 2-wave

| $Y_0/G$ | $v_p^*$ split (cm/s) | $v_p^*$ original (cm/s) | Shift |
|---------|---------------------|------------------------|-------|
| 0.005 | 86 077 | 86 273 | +0.23% |
| 0.055 | 83 658 | 85 809 | +2.51% |
| 0.105 | 81 219 | 85 305 | +4.79% |
| 0.155 | 78 760 | 84 759 | +7.08% |
| 0.205 | 76 281 | 84 170 | +9.37% |

---

## 6. Phase 5 -- Physical Interpretation of Regimes

### Overdriven regime

When $v_\mathrm{piston} \geq U_{se}$, the piston velocity equals or exceeds
the elastic precursor speed. The two-wave structure is impossible because the
piston would "outrun" the elastic wave. In practice, a single overdriven shock
forms, but this is outside the scope of the two-wave model. With the
Mie-Gruneisen EOS, $U_{se}$ is so high that overdriven conditions fall outside
the plotted range.

### No-root regime

When $f_\mathrm{shock}$ has no sign change in $(v_\mathrm{piston}, U_{se})$,
no self-consistent two-wave shock exists. This can happen because:

- At high $v_p$, the Rankine-Hugoniot momentum pressure $P_2(U_s)$ is
  dominated by the EOS pressure $P_\mathrm{eos}(\rho_2, e_2)$, and the
  equation $P_2 = P_\mathrm{eos}$ has no solution.
- Near the MG Hugoniot singularity ($\mu \to 1/s$), the EOS pressure
  diverges, making $f_\mathrm{shock} < 0$ everywhere.

### Mie-Gruneisen Hugoniot singularity

The Hugoniot reference pressure is:

$$P_H(\mu) = \frac{\rho_0\,C_0^2\,\mu}{(1 - s\mu)^2}$$

This diverges at $\mu = 1/s$ (i.e., $\rho = \rho_0\,s/(s-1)$). For our
parameters ($s = 1.34$), the singularity is at $\mu = 0.746$, corresponding
to $\rho = 10.98$ g/cm$^3$ ($3.94\times$ reference density).

This singularity:

- Sets a **maximum density** beyond which the MG EOS is undefined.
- Creates the **no-root region** at extreme piston velocities.
- Makes $f_\mathrm{shock}$ develop NaN/Inf values near the singularity,
  preventing root-finding.
- Is a consequence of the linear $U_s = C_0 + s\,u_p$ Hugoniot fit,
  which is only an empirical approximation valid for moderate compressions.

---

## 7. Phase 6 -- Ideal Gas EOS Regime Maps

**Script**: [`scripts/regime_map.py`](scripts/regime_map.py) (extended with
`_eos_ig` and `eos_fn` parameter threading)

To isolate the effect of the MG Hugoniot singularity, we replaced the EOS with
an ideal gas:

$$P = (\gamma - 1)\,\rho\,e, \qquad \gamma = \Gamma_0 + 1 = 3$$

This preserves the thermal-pressure coefficient $\Gamma_0\,\rho$ but eliminates
the Hugoniot reference curve ($P_H = 0$, $e_H = 0$).

### IG energy-split mode

![IG energy-split regime map](assets/regime_map_ig_split.png)

### IG original mode

![IG original mode regime map](assets/regime_map_ig_orig.png)

### MG vs IG boundary comparison

![MG vs IG boundaries](assets/regime_map_ig_vs_mg.png)

### Regime statistics comparison

| Regime | MG split (%) | MG orig (%) | IG split (%) | IG orig (%) |
|--------|-------------|------------|-------------|-------------|
| No plastic shock | 34.4 | 35.1 | 22.3 | 23.7 |
| Physical 2-wave | 11.7 | 12.5 | 35.3 | 35.9 |
| Strong shock only | 46.2 | 45.2 | **0.0** | **0.0** |
| No root | 7.7 | 7.1 | 38.5 | 38.3 |
| Overdriven | 0.0 | 0.0 | 3.9 | 2.1 |

### Key findings

1. **The "strong shock only" region vanishes entirely** with the ideal gas EOS.
   Without $P_H(\mu)$, $f_\mathrm{shock}$ does not develop the intermediate root
   structure that produces the unphysical strong-shock solution. The transition
   goes directly from physical 2-wave to no-root.

2. **The "no root" region expands, not shrinks.** This was the key surprise.
   Without the cold pressure $P_H$, the EOS thermal pressure $(\gamma-1)\rho e$
   exceeds the Rankine-Hugoniot momentum pressure $P_2$ across the entire $U_s$
   range at moderate-to-high piston velocities. The function is well-behaved (no
   NaN/Inf) but simply never crosses zero. The two-wave structure genuinely
   cannot satisfy both the RH jump conditions and the ideal gas EOS at those
   velocities.

3. **The physical 2-wave wedge roughly doubles** ($v_p^*/C_0 \approx 0.3$
   instead of $\approx 0.16$), because the ideal gas EOS lacks the stiff
   cold-curve pressure that flips $f_\mathrm{shock}$ negative near $U_{se}$.

4. **The "Overdriven" region becomes visible** in the bottom-right corner, where
   $v_p$ exceeds $U_{se}$. With MG, the elastic precursor speed is much higher
   (due to $P_H$), pushing this region beyond the plotted range.

5. **The energy-split effect persists** with the ideal gas EOS, because elastic
   energy subtraction is EOS-independent.

---

## 8. Phase 7 -- Wave Speed Ratio Analysis

**Script**: [`scripts/regime_map.py`](scripts/regime_map.py)
(functions `solve_Us`, `build_ratio_grid`, `plot_ratio_map`)

Inside the physical 2-wave wedge, we computed

$$\log_{10}\!\left(\frac{U_{se}}{U_s} - 1\right)$$

which measures the fractional speed excess of the elastic precursor over the
plastic shock. On this scale:
- $-3$ means the two waves differ by 0.1%
- $-1$ means 10% difference
- $0$ means 100% (factor of 2)

### 2x2 panel (all four EOS/mode combinations)

![Ratio panel](assets/ratio_panel.png)

### Individual maps

| | Energy-split | Original |
|---|---|---|
| **MG** | ![MG split ratio](assets/ratio_mg_split.png) | ![MG orig ratio](assets/ratio_mg_orig.png) |
| **IG** | ![IG split ratio](assets/ratio_ig_split.png) | ![IG orig ratio](assets/ratio_ig_orig.png) |

### Observations

- **Mie-Gruneisen**: the separation is small, $\log_{10}(U_{se}/U_s - 1)$
  ranges from about $-3$ to $-1$ (0.1% to 10%). The plastic shock speed is
  always close to the elastic precursor speed.

- **Ideal gas**: the separation is much larger, reaching
  $\log_{10}(U_{se}/U_s - 1) \approx 1$ (factor of 10 difference). At higher
  $v_\mathrm{piston}$ the plastic shock slows dramatically relative to the
  elastic precursor. This reflects the fact that, without the cold-curve
  stiffness, the plastic wave does not "keep up" with the elastic one.

- In both EOS cases, the ratio increases smoothly from left to right (higher
  $v_\mathrm{piston}$), indicating greater wave separation at stronger loading.

---

## 9. Phase 8 -- Ideal Gas Profile Comparison

**Script**: [`scripts/ig_mode_comparison.py`](scripts/ig_mode_comparison.py)

To visualize the actual differences between original and energy-split modes, we
compared density, stress, and velocity profiles for three ideal-gas cases
chosen to maximize the visible difference:

| Case | $Y_0/G$ | $v_\mathrm{piston}$ (cm/s) | $t$ ($\mu$s) |
|------|---------|---------------------------|--------------|
| 1 | 0.20 | 40 000 | 1.5 |
| 2 | 0.30 | 80 000 | 0.8 |
| 3 | 0.50 | 120 000 | 0.5 |

### Case 1: $Y_0/G = 0.20$, $v_p = 40\,000$ cm/s

![Case 1 comparison](assets/ig_case1_comparison.png)

### Case 2: $Y_0/G = 0.30$, $v_p = 80\,000$ cm/s

![Case 2 comparison](assets/ig_case2_comparison.png)

### Case 3: $Y_0/G = 0.50$, $v_p = 120\,000$ cm/s

![Case 3 comparison](assets/ig_case3_comparison.png)

### Combined 3x3 panel

![IG comparison panel](assets/ig_comparison_panel.png)

### Observations

- At $Y_0/G = 0.50$, $v_p = 120\,000$ cm/s (Case 3), the relative difference
  in $P_2$ between modes is ~36%, and $U_s$ differs by ~9%.
- The energy-split mode consistently produces **higher compression** (higher
  $\rho_2$, more negative $\sigma_x$, lower $U_s$) because subtracting elastic
  energy from the EOS input lowers the thermal pressure, requiring more
  compression to satisfy the momentum jump condition.
- The elastic precursor state ($\rho_Y$, $P_Y$, $v_Y$) also differs between
  modes, contributing to distinct precursor positions.
- These differences grow with $Y_0/G$ because the elastic energy fraction
  $e_\mathrm{el}/e_\mathrm{total}$ increases with yield strength.

---

## 10. Phase 9 -- MG-to-Ideal-Gas Parameter Mapping

The Mie-Gruneisen EOS reduces to an ideal gas when the Hugoniot reference curve
is eliminated. Setting $C_0 = 0$ gives:

$$P_H = 0, \quad e_H = 0$$

so:

$$P = P_H + \Gamma_0\,\rho\,(e - e_H) = \Gamma_0\,\rho\,e = (\gamma - 1)\,\rho\,e$$

with $\gamma = \Gamma_0 + 1$. The $s$ parameter becomes irrelevant since
$P_H = 0$.

### Parameters for MG behaving as ideal gas

```python
rho_0   = 2.79          # g/cm^3       reference density
C_0     = 0.0           # cm/s         kills Hugoniot -> ideal gas
s       = 1.34          #              irrelevant (P_H = 0)
Gamma_0 = 2.0           #              gives gamma = Gamma_0 + 1 = 3
G       = 2.86e11       # dyn/cm^2     shear modulus
Y_0     = 1.43e11       # dyn/cm^2     yield stress  (Y0/G = 0.50)

v_piston  = 120000.0    # cm/s         piston velocity
e_initial = 0.0         # erg/g        initial internal energy

t     = 0.5e-6          # s            snapshot time
x_max = 0.29            # cm           domain length
```

This provides a way to test ideal-gas-like behavior using the existing MG solver
without modifying the EOS code.

---

## 11. Conclusions

1. **The loss of the weak-shock root is tied to the Hugoniot reference curve**
   and its singularity at $\mu = 1/s$. Both modes lose the physical solution
   above $v_p/C_0 \approx 0.163$.

2. **Energy splitting shifts the critical boundary to lower
   $v_\mathrm{piston}$**. At $Y_0/G = 0.15$, the shift is large enough that
   $v_p = 80\,000$ cm/s falls outside the physical 2-wave region in
   energy-split mode but remains inside it in original mode.

3. **The physical 2-wave window is a triangular wedge** bounded below by
   $v_Y(Y_0/G)$ (diagonal line) and above by $v_p^*(Y_0/G)$ (near-horizontal
   curve). As $Y_0/G$ increases, $v_Y$ rises and $v_p^*$ drops, squeezing the
   wedge closed.

4. **At extreme velocities ($v_p/C_0 \gtrsim 0.65$), $f_\mathrm{shock}$ loses
   all roots**, creating the "no root" region. This is caused by the MG
   Hugoniot singularity.

5. **The ideal gas EOS eliminates the Hugoniot singularity.** With
   $P = (\gamma - 1)\rho e$, there is no $(1 - s\mu)^2$ denominator.

6. **The "strong shock only" region vanishes entirely** with the ideal gas EOS.
   The transition goes directly from physical 2-wave to no-root, confirming
   that the strong-shock-only artifact is specific to the MG cold curve.

7. **The "no root" region expands with ideal gas**, replacing the MG
   "strong shock only" failure. Neither EOS produces valid two-wave solutions
   above $v_p/C_0 \approx 0.3$ (IG) or $\approx 0.16$ (MG). The loss of
   *any* solution at high $v_p$ is a structural property of the two-wave
   elastoplastic model.

8. **The physical 2-wave wedge roughly doubles** under ideal gas EOS,
   extending to $v_p/C_0 \approx 0.3$ instead of $\approx 0.16$.

9. **Wave speed separation differs dramatically**: MG gives 0.01--10% ratio
   ($U_{se}/U_s - 1$) while ideal gas gives up to a factor of 10.

10. **The energy-split effect is EOS-independent**, always shifting the 2-wave
    boundary to lower $v_\mathrm{piston}$. The solver fix (scanning all roots,
    picking the weakest, warning on strong-shock-only) is the correct approach
    for both EOS types.

### Practical guidance

When using the energy-split formulation, ensure:

$$v_Y(Y_0/G) < v_\mathrm{piston} < v_p^*(Y_0/G)$$

The $v_p^*$ boundary can be pre-computed for a given material. The solver issues
`warnings.warn()` when only the strong-shock root is found, alerting users that
the result may be unphysical.

---

## Scripts Reference

| Script | Purpose |
|--------|---------|
| [`diagnose_energy_split.py`](scripts/diagnose_energy_split.py) | Scans $f_\mathrm{shock}(U_s)$, finds all roots, checks physical constraints, generates diagnostic plots |
| [`regime_map.py`](scripts/regime_map.py) | Generates all regime maps (MG and IG), boundary curves, ratio maps, and the detailed report |
| [`ig_mode_comparison.py`](scripts/ig_mode_comparison.py) | Compares density/stress/velocity profiles between original and energy-split modes for ideal gas cases |
| [`elastoplastic_piston_solver.py`](scripts/elastoplastic_piston_solver.py) | Modified solver with robust multi-root scanning in `_find_physical_shock_speed` |
