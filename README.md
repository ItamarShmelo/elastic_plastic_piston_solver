# Elastoplastic Piston Test Problem

The problem is a one dimensional hydrodynamic elastoplastic test
problem. It consists of a 1D medium between $\left[0,L\right]$, starting
at rest i.e velocity 0 and specific internal energy $e=e_{initial}$,
with a piston at $x=0$ moving the in the direction of the positive
$x$-axis with velocity $U_{piston}$. The material is a perfectly
elastic-plastic material i.e. constant shear modulus $G$ and yield
strength $Y_{0}$ it has a Mie-Gruneisen equation of state (see
[Mie-Gruneisen EOS](#mie-gruneisen-eos)), or an ideal-gas EOS (see
[Ideal-Gas EOS](#ideal-gas-eos)).

Depending on the input parameters, the solver returns one of three
wave structures:

| Wave structure | Condition | Description |
|---|---|---|
| **Elastic wave** | $v_{piston} \leq v^{Y}$ | Piston velocity is at or below the elastic-limit particle velocity; only an elastic precursor propagates. A full 0-to-E Rankine-Hugoniot solve is performed with the actual sub-yield deviatoric stress $S_x^E = \tfrac{4}{3}G\ln(\rho_0/\rho_E)$ and work-integral elastic energy. |
| **Plastic shock and elastic wave** | $U_{s,02} \leq U_{se}$ | Standard two-wave regime: a plastic shock trailed by an elastic precursor. The direct 0-to-2 shock speed does not overrun the elastic precursor. |
| **Plastic wave** | $U_{s,02} > U_{se}$ | A single plastic shock with no separate elastic precursor. The direct 0-to-2 shock speed overruns the elastic precursor, so the Rankine-Hugoniot jump is taken directly from the initial state to the shocked state. |

## Example: Aluminum Piston

The following example solves the elastoplastic piston problem for
**Aluminum** (CGS units) with a piston velocity of 5000 cm/s, evaluated
at $t = 2\;\mu s$.

```python
import numpy as np
from elastoplastic_piston_solver import ElastoplasticPistonSolver

# Material parameters (Aluminum, CGS)
solver = ElastoplasticPistonSolver(
    rho_0=2.79,        # g/cm^3    — reference density
    C_0=5.33e5,        # cm/s      — bulk sound speed
    s=1.34,            #           — Hugoniot slope coefficient
    Gamma_0=2.0,       #           — Gruneisen parameter
    G=2.86e11,         # dyn/cm^2  — shear modulus
    Y_0=2.6e9,         # dyn/cm^2  — yield strength
    e_initial=0.0,     # erg/g     — initial specific internal energy
    v_piston=5000.0,   # cm/s      — piston velocity
)

# Evaluate at t = 2 microseconds
t = 2.0e-6  # seconds
x_max = 1.2 * solver.U_se * t
x = np.linspace(0.0, x_max, 1000)

result = solver.solve(t, x)
```

The returned `result` dictionary contains arrays for `"density"`,
`"pressure"`, `"velocity"`, `"energy"`, `"Sx"` (deviatoric stress), and
`"stress"` (total axial stress $\sigma_x = S_x - P$), as well as the
scalar positions `"shock_location"` and `"elastic_precursor_location"`,
and the `"wave_structure"` enum value (`WaveStructure`).

When a wave front is absent for the given regime, its location is
`np.nan`:

| Wave structure | `shock_location` | `elastic_precursor_location` |
|---|---|---|
| `ELASTIC_WAVE` | `np.nan` | finite |
| `PLASTIC_SHOCK_AND_ELASTIC_WAVE` | finite | finite |
| `PLASTIC_WAVE` | finite | `np.nan` |

See the [Example Figures](#example-figures) section at the end of this
document for stress and velocity profile plots across all regimes and
EOS types.  To reproduce the aluminum plots, run:

```bash
python examples/example_aluminum_piston.py
```

## Energy-Split Mode

The solver supports an optional **energy-split mode** that decomposes the
internal energy into thermal and elastic contributions:

$$e = e_{th} + e_{el}$$

The elastic energy is computed from the **work integral** of the elastic
deviatoric power $\boldsymbol{S}:\boldsymbol{D}^{el}$:

$$e_{el}(\rho) = \frac{4G}{3\rho_0}\left[1 - (1+a)\,e^{-a}\right], \qquad a = \ln\frac{\rho}{\rho_0}$$

This is the integral $\int_{\rho_0}^{\rho} (-S_x)\,d\rho'/\rho'^2$ with
$S_x = \tfrac{4}{3}G\ln(\rho_0/\rho')$.  It matches the thermal energy
PDE used in hydrocodes ($\rho\,de_{th}/dt = \rho\mathcal{P}_{plast} - P\nabla\cdot\vec{u}$),
where only pressure work and plastic dissipation heat the material.

During plastic deformation the deviatoric stress is constant at yield
($\boldsymbol{D}^{el} = 0$), so $e_{el}$ is frozen at its yield-point
value.

> **Note:** The previous formulation used the state function
> $e_{el} = \boldsymbol{S}:\boldsymbol{S}/(4\rho G)$, which differs from
> the work integral for finite elastic strains ($Y_0/(2G) \gtrsim 0.05$).
> For typical metals ($Y_0/(2G) \sim 10^{-3}$) the two are equivalent to
> leading order.

In this mode the EOS receives only the thermal part, $P(e_{th}, \rho)$,
rather than the total energy.  Enable it by passing `energy_split=True`:

```python
solver = ElastoplasticPistonSolver(
    rho_0=2.79, C_0=5.33e5, s=1.34, Gamma_0=2.0,
    G=2.86e11, Y_0=2.6e9,
    e_initial=0.0, v_piston=5000.0,
    energy_split=True,   # <-- enable energy-split mode
)
```

When `energy_split=True` the returned dictionary from `solver.solve()`
includes two additional arrays, `"e_thermal"` and `"e_elastic"`, alongside
the usual fields.  The `"energy"` key continues to hold the total energy.

### Modified equations

Let $e_{el}^Y = \frac{4G}{3\rho_0}\left[1 - \left(1 + \frac{Y_0}{2G}\right)e^{-Y_0/(2G)}\right]$ denote the work-integral elastic energy at yield.

The yield-energy root-finding equation becomes:

$$f_{Y}(e_{th}^{Y}) = e_{th}^{Y} + e_{el}^Y - e_0 - \frac{1}{2\rho^{Y}\rho_0}\left(P^{Y}(e_{th}^{Y},\rho^{Y}) + \tfrac{2}{3}Y_0\right)(\rho^{Y} - \rho_0)$$

The shocked thermal energy uses $e_{el}$ frozen at yield (since $\boldsymbol{D}^{el} = 0$ during plastic deformation):

$$e_{th,2} = e_{th}^{Y} + e_{el}^Y + \frac{1}{2\rho^{Y}\rho_2}\left(P^{Y} + P_2 + \tfrac{4}{3}Y_0\right)(\rho_2 - \rho^{Y}) - e_{el}^Y$$

which simplifies to:

$$e_{th,2} = e_{th}^{Y} + \frac{1}{2\rho^{Y}\rho_2}\left(P^{Y} + P_2 + \tfrac{4}{3}Y_0\right)(\rho_2 - \rho^{Y})$$

See the [Derivation of the Analytic Solution](#derivation-of-the-analytic-solution)
section below for the full derivation.

### Mode comparison

For typical metals ($Y_0/G \sim 10^{-3}$) the two modes give nearly
identical results ($< 0.1\%$ relative difference).  The effect becomes
significant when $Y_0/G \gtrsim 0.05$.

A detailed quantitative comparison for Aluminum, Copper, and a
high-strength test case is available in
[`examples/energy_split_comparison.md`](examples/energy_split_comparison.md).

To regenerate the comparison plots and report, run:

```bash
python examples/example_mode_comparison.py
```

See the [Example Figures](#example-figures) section for the
high-strength comparison plot.

## Derivation of the Analytic Solution

We start the deviatoric stress, below the yield limit.

$$\frac{dS_{x}}{dt}=\frac{4}{3}G\frac{\partial v}{\partial x}$$

Where

$$S_{x}=S_{xx},\ S_{y}\equiv S_{yy},\ S_{z}\equiv S_{zz}$$

Which when integrated for uni-axial motion gives (since
$\frac{1}{\rho}\frac{d\rho}{dt}=-\frac{\partial v}{\partial x}$)

$$S_{x}=\frac{4}{3}G\ln\frac{\rho_{0}}{\rho}$$

This holds until
$\boldsymbol{S}:\boldsymbol{S}\leq\frac{2}{3}Y_{0}^{2}$ or

$$S_{x}^{2}+S_{y}^{2}+S_{z}^{2}\leq\frac{2}{3}Y_{0}^{2}$$

since there are no shear stresses. Using $S_{z}=S_{y}=-\frac{1}{2}S_{x}$,
since $S_{x}+S_{y}+S_{z}=0$ and symmetry. We get that at and above the yield
point

$$\left(S_{x}^{Y}\right)^{2}+\frac{1}{4}\left(S_{x}^{Y}\right)^{2}+\frac{1}{4}\left(S_{x}^{Y}\right)^{2}\leq\frac{2}{3}Y_{0}^{2}$$

Setting $S=S_{x}$

$$\frac{3}{2}\left(S^{Y}\right)^{2}\leq\frac{2}{3}Y_{0}^{2}$$

$$\left(S^{Y}\right)^{2}\leq\frac{4}{9}Y_{0}^{2}$$

$$S^{Y}=-\frac{2}{3}Y_{0}$$

Thus after reaching the yield point

$$S^{Y}=-\frac{2}{3}Y_{0}=\frac{4}{3}G\ln\frac{\rho_{0}}{\rho^{Y}}$$

$$\rho^{Y}=\rho_{0}e^{\frac{Y_{0}}{2G}}$$

The Hugoniot relations (from Wilkins).

$$U_{shock}^{2}=-\frac{1}{\rho_{0}}\frac{\sigma_{1}-\sigma_{0}}{\left(1-\frac{\rho_{0}}{\rho_{1}}\right)}=\frac{\rho_{1}\left(\sigma_{0}-\sigma_{1}\right)}{\rho_{0}\left(\rho_{1}-\rho_{0}\right)}$$

where, $\sigma\equiv\sigma_{x}$ and $\sigma=S-P$ since
$\sigma_{0}\approx0$,

$$\sigma_{1}=\sigma^{Y}=-\frac{2}{3}Y_{0}-P^{Y}=-\left(\frac{2}{3}Y_{0}+P^{Y}\right)$$

and $\rho^{Y}$ we get that

$$U_{se}^{2}=\frac{\rho^{Y}\sigma^{Y}}{\rho_{0}(\rho_{0}-\rho^{Y})}=-\frac{\rho^{Y}\left(P^{Y}+\frac{2}{3}Y_{0}\right)}{\rho_{0}(\rho_{0}-\rho^{Y})}$$

The internal energy from the Hugoniot relation

$$e_{1}-e_{0}=-\frac{1}{2\rho_0}\left(\sigma_{1}+\sigma_{0}\right)\left(1-\frac{\rho_{0}}{\rho_{1}}\right)$$

Or

$$e^{Y}=e_{0}+\frac{1}{2\rho^{Y}\rho_{0}}\left(P^{Y}\left(e^{Y},\rho^{Y}\right)+\frac{2}{3}Y_{0}\right)\left(\rho^{Y}-\rho_{0}\right)$$

This is an implicit equation for $e^{Y}$ since
$P^{Y}\left(e^{Y},\rho^{Y}\right)$ via the Mie-Gruneisen EOS. After
solving the above we can get the particle velocity (from the Hugoniot
relation

$$v^{Y}=\frac{\rho^{Y}-\rho_{0}}{\rho^{Y}}U_{se}$$

To solve for the plastic shock we solve the Hugoniot system

$$P_{2}=P^{Y}+\rho^{Y}\left(U_{s}-v^{Y}\right)\left(v_{2}-v^{Y}\right)$$

$$\rho_{2}=\rho^{Y}\frac{U_{s}-v^{Y}}{U_{s}-v_{2}}$$

$$e_{2}=e^{Y}+\frac{1}{2\rho^{Y}\rho_{2}}\left(P^{Y}+P_{2}+\frac{4}{3}Y_{0}\right)\left(\rho_{2}-\rho^{Y}\right)$$

since $\sigma_{2}=-P_{2}-\frac{2}{3}Y_{0}$ since we are already past the
yield point at this time.

$$v_{2}=v_{piston}$$

### One-wave regime: direct 0-to-2 jump

When the direct 0-to-2 shock speed $U_{s,02}$ exceeds the elastic
precursor speed $U_{se}$, the plastic shock overruns the precursor and
the two-wave structure is not viable. The solver uses a single plastic
wave, computing $U_s$ and the shocked state via the Rankine-Hugoniot
jump directly from the initial state 0 to state 2.

State 0 has zero deviatoric stress ($S_x^0 = 0$, $\sigma_0 \approx 0$),
while state 2 is at yield ($S_x^2 = -\tfrac{2}{3}Y_0$). The jump
conditions are:

$$\rho_2 = \frac{\rho_0 U_s}{U_s - v_{piston}}$$

$$P_2 = \rho_0 U_s v_{piston} - \frac{2}{3}Y_0$$

$$e_2 = e_0 + \frac{1}{2\rho_0\rho_2}\left(P_2 + \frac{2}{3}Y_0\right)\left(\rho_2 - \rho_0\right)$$

The residual for root-finding is:

$$f_S(U_s) = P_2(U_s) - P_{eos}\left(e_2(U_s),\;\rho_2(U_s)\right) = 0$$

When `energy_split=True`, the EOS receives
$e_{th,2} = e_2 - e_{el}^Y$ (elastic energy frozen at yield).

### Elastic wave regime: 0-to-E jump

When $v_{piston} \leq v^Y$, the material does not reach yield. The
elastic wave is a shock from state 0 (undisturbed) to state E
(sub-yield compression). Unlike the two-wave regime, the deviatoric
stress is **not** at yield:

$$S_x^E = \frac{4}{3}G\ln\frac{\rho_0}{\rho_E}$$

The elastic energy is computed from the work integral at the sub-yield
density $\rho_E$:

$$e_{el}^E = \frac{4G}{3\rho_0}\left[1 - (1+a_E)\,e^{-a_E}\right], \qquad a_E = \ln\frac{\rho_E}{\rho_0}$$

At yield ($\rho_E = \rho^Y$, $a_E = Y_0/(2G)$), this gives the elastic
energy at the elastic limit, ensuring continuity between the elastic and
two-wave regimes.

Given elastic shock speed $U_{se}$, the Rankine-Hugoniot jump conditions
from state 0 to state E are:

$$\rho_E = \frac{\rho_0 U_{se}}{U_{se} - v_{piston}}$$

$$P_E = \rho_0 U_{se} v_{piston} + S_x^E$$

$$e_E = e_0 + \frac{1}{2\rho_0\rho_E}\left(P_E - S_x^E\right)\left(\rho_E - \rho_0\right)$$

The residual for root-finding is:

$$f_{E}(U_{se}) = P_E - P_{eos}\left(e_{eos},\;\rho_E\right) = 0$$

where $e_{eos} = e_E$ when `energy_split=False`, or
$e_{eos} = e_E - e_{el}^E$ when `energy_split=True`.

The scan range is $v_{piston} \cdot 1.001$ to
$\max(2 U_{se}^{yield},\; 10\, v_{piston})$, with the **largest** root
selected (weakest compression, physically correct solution).

### Equivalent plastic strain

The equivalent (von Mises) plastic strain is defined by

$$\bar\varepsilon^{pl} = \int_0^t \sqrt{\tfrac{2}{3}\,\boldsymbol{D}^{pl}:\boldsymbol{D}^{pl}}\;dt$$

The hypoelastic J2 plasticity model evolves the deviatoric stress via

$$\overset{\nabla}{\boldsymbol{S}} = 2G\left(\bar{\boldsymbol{D}} - \boldsymbol{D}^{pl}\right)$$

where $\bar{\boldsymbol{D}} = \boldsymbol{D} - \tfrac13(\nabla\!\cdot\!\boldsymbol{v})\,\boldsymbol{I}$
is the deviatoric part of the strain-rate.  The plastic strain-rate
$\boldsymbol{D}^{pl}$ is therefore a deviatoric tensor, and the
equivalent plastic strain must be computed from this deviatoric
quantity, not from the total post-yield axial strain.

For uniaxial strain along $y$, the total strain-rate is
$\boldsymbol{D} = \mathrm{diag}(0,\;\partial v/\partial y,\;0)$, so

$$\bar{\boldsymbol{D}}
= \mathrm{diag}\!\left(-\tfrac13\frac{\partial v}{\partial y},\;
\tfrac23\frac{\partial v}{\partial y},\;
-\tfrac13\frac{\partial v}{\partial y}\right)$$

In the perfectly plastic region the deviatoric stress is clamped at
yield ($\overset{\nabla}{\boldsymbol{S}}=0$, no spin in 1-D), so
$\boldsymbol{D}^{pl} = \bar{\boldsymbol{D}}$ and in particular

$$D^{pl}_{yy} = \tfrac23\frac{\partial v}{\partial y}$$

The double contraction is
$\boldsymbol{D}^{pl}:\boldsymbol{D}^{pl}
= \tfrac23\left(\frac{\partial v}{\partial y}\right)^2$,
giving

$$\dot{\bar\varepsilon}^{pl}
= \sqrt{\tfrac{2}{3}\cdot\tfrac{2}{3}}\left|\frac{\partial v}{\partial y}\right|
= \tfrac{2}{3}\left|\frac{\partial v}{\partial y}\right|$$

Using the continuity equation
$\partial v/\partial y = -(\rho)^{-1}\,d\rho/dt$, for monotone
compression ($d\rho/dt > 0$):

$$\dot{\bar\varepsilon}^{pl}
= \frac{2}{3}\frac{1}{\rho}\frac{d\rho}{dt}$$

Integrating from $\rho^Y$ to $\rho$ gives the equivalent plastic
strain as a function of density:

$$\bar\varepsilon^{pl}(\rho)=
\begin{cases}
0, & \rho \le \rho^Y \\[4pt]
\dfrac{2}{3}\ln\!\left(\dfrac{\rho}{\rho^Y}\right), & \rho > \rho^Y
\end{cases}$$

In the shocked state 2:

$$\bar\varepsilon^{pl}_2 = \frac{2}{3}\ln\!\left(\frac{\rho_2}{\rho^Y}\right)$$

## Solver Algorithm

1.  Define problem parameters

$$\text{EOS} :\rho_{0},C_{0},s,\Gamma_{0}$$

$$\text{Elastoplastic} :G,Y_{0}$$

$$\text{Initial Conditions} :e_{initial},v_{piston}$$

2.  Calculate $\rho^{Y}$ via

$$\rho^{Y}=\rho_{0}e^{\frac{Y_{0}}{2G}}$$

3.  Calculate $e^{Y},p^{Y}$ by finding the root of the function

$$f_{Y}\left(e^{Y}\right)=e^{Y}-e_{0}-\frac{1}{2\rho^{Y}\rho_{0}}\left(P^{Y}\left(e^{Y},\rho^{Y}\right)+\frac{2}{3}Y_{0}\right)\left(\rho^{Y}-\rho_{0}\right)$$

or with the energy-split mode, solve for $e_{th}^{Y}$:

$$f_{Y}\left(e_{th}^{Y}\right)=e_{th}^{Y}+e_{el}^Y-e_{0}-\frac{1}{2\rho^{Y}\rho_{0}}\left(P^{Y}\left(e_{th}^{Y},\rho^{Y}\right)+\frac{2}{3}Y_{0}\right)\left(\rho^{Y}-\rho_{0}\right)$$

where $e_{el}^Y = \frac{4G}{3\rho_0}\left[1 - \left(1 + \frac{Y_0}{2G}\right)e^{-Y_0/(2G)}\right]$.

4.  Calculate $U_{se}$

$$U_{se}^{2}=-\frac{\rho^{Y}\left(P^{Y}+\frac{2}{3}Y_{0}\right)}{\rho_{0}(\rho_{0}-\rho^{Y})}$$

5.  Calculate $v^{Y}$ via

$$v^{Y}=\frac{\rho^{Y}-\rho_{0}}{\rho^{Y}}U_{se}$$

6.  **Regime classification.**

    - If $v_{piston} \leq v^Y$: **elastic wave** — perform a full 0-to-E
      shock solve (see [Elastic wave regime](#elastic-wave-regime-0-to-e-jump))
      with the actual sub-yield $S_x^E$ and $e_{el}^E$, then skip steps 7–8.
    - Otherwise, solve the 0-to-2 Rankine-Hugoniot residual to obtain the
      direct shock speed $U_{s,02}$ (scan from $v_{piston} \cdot 1.001$
      to $\max(2 U_{se},\; 5\, v_{piston})$, pick the largest root).
    - If $U_{s,02} > U_{se}$: **plastic wave** (one-wave regime)
      — continue to step 7b using $U_s = U_{s,02}$.
    - If $U_{s,02} \leq U_{se}$: **plastic shock and elastic wave**
      (two-wave regime) — scan the Y-to-2 residual below $U_{se}$
      and continue to step 7a.

7a. **Two-wave U_s solve (Y-to-2 jump).** Pick the largest root below
    $U_{se}$ (the weak-shock root nearest the elastic precursor). Then
    compute state 2 from the Y-to-2 jump conditions:

$$P_{2}=P^{Y}+\rho^{Y}\left(U_{s}-v^{Y}\right)\left(v_{2}-v^{Y}\right)$$

$$\rho_{2}=\rho^{Y}\frac{U_{s}-v^{Y}}{U_{s}-v_{2}}$$

$$e_{2}=e^{Y}+\frac{1}{2\rho^{Y}\rho_{2}}\left(P^{Y}+P_{2}+\frac{4}{3}Y_{0}\right)\left(\rho_{2}-\rho^{Y}\right)$$

or with the energy-split mode, use $e_{th,2}$ in the EOS (elastic energy frozen at yield):

$$e_{th,2}=e_{th}^{Y}+\frac{1}{2\rho^{Y}\rho_{2}}\left(P^{Y}+P_{2}+\frac{4}{3}Y_{0}\right)\left(\rho_{2}-\rho^{Y}\right)$$

7b. **One-wave U_s solve (0-to-2 jump).** Build the direct 0-to-2
    residual and scan from $v_{piston} \cdot 1.001$ to
    $\max(2 U_{se},\; 5 v_{piston})$. Pick the largest root. Then
    compute state 2 from the 0-to-2 jump conditions:

$$\rho_{2}=\frac{\rho_0 U_s}{U_s - v_{piston}},\quad P_2 = \rho_0 U_s v_{piston} - \tfrac{2}{3}Y_0$$

$$e_2 = e_0 + \frac{1}{2\rho_0\rho_2}\left(P_2 + \tfrac{2}{3}Y_0\right)\left(\rho_2 - \rho_0\right)$$

8.  **Piecewise-constant profile.** When given a grid
    $x_1, x_2, \ldots$ and a time $t$:

For the **plastic shock and elastic wave** (two-wave) regime:

$$
X(x_i) = \left\lbrace \begin{array}{ll}
X_2 & x_i < U_s t \\
X^Y & U_s t \leq x_i \leq U_{se} t \\
X_{\text{initial}} & U_{se} t < x_i
\end{array} \right.
$$

For the **elastic wave** regime, only the elastic precursor propagates
(no shocked region). For the **plastic wave** regime, only the plastic
shock propagates (no elastic precursor region).

## Equation of State

The solver supports two EOS configurations, selected through the
constructor.  Exactly one must be provided; mixing or omitting both
raises `ValueError`.

### Mie-Gruneisen EOS

Provide `C_0`, `s`, and `Gamma_0` (all three required together):

```python
solver = ElastoplasticPistonSolver(
    rho_0=2.79, C_0=5.33e5, s=1.34, Gamma_0=2.0,
    G=2.86e11, Y_0=2.6e9, e_initial=0.0, v_piston=5000.0,
)
```

The Mie-Gruneisen equation of state:

$$P\left(e,\rho\right)=P_{H}\left(\rho\right)+\Gamma_{0}\rho\left(e-e_{H}\left(\rho\right)\right)$$

for compression $\rho\geq\rho_{0}$

$$P_{H}\left(\rho\right)=\frac{\rho_{0}C_{0}^{2}\mu}{\left(1-s\mu\right)^{2}}$$

where

$$\mu=1-\frac{\rho_{0}}{\rho}$$

is the compression and

$$e_{H}\left(\rho\right)=\frac{1}{2}P_{H}\left(\rho\right)\frac{\mu}{\rho_{0}}=\frac{1}{2}C_{0}^{2}\frac{\mu^{2}}{\left(1-s\mu\right)^{2}}$$

So together

$$P\left(e,\rho\right)=\frac{\rho_{0}C_{0}^{2}\mu}{\left(1-s\mu\right)^{2}}+\Gamma_{0}\rho\left(e-\frac{C_{0}^{2}\mu^{2}}{2\left(1-s\mu\right)^{2}}\right)$$

<div align="center">

| Parameter    | Description                              |
|--------------|------------------------------------------|
| $\rho_{0}$   | Reference density                        |
| $C_{0}$      | Bulk sound speed at reference state      |
| $s$          | Hugoniot slope coefficient               |
| $\Gamma_{0}$ | Gruneisen parameter (assumed constant)   |

</div>

### Ideal-Gas EOS

Provide `gamma_ideal_gas` only (mutually exclusive with `C_0`/`s`/`Gamma_0`):

```python
solver = ElastoplasticPistonSolver(
    rho_0=2.79, gamma_ideal_gas=3.0,
    G=2.86e11, Y_0=2.6e9, e_initial=0.0, v_piston=5000.0,
)
```

The ideal-gas EOS $P = (\gamma - 1)\rho e$ is implemented internally by
mapping to Mie-Gruneisen constants that zero out the Hugoniot reference
curve:

$$C_0 = 0, \quad s = 0, \quad \Gamma_0 = \gamma - 1$$

With $C_0 = 0$, both $P_H$ and $e_H$ vanish identically, and the
Mie-Gruneisen EOS reduces to $P = \Gamma_0 \rho\, e = (\gamma-1)\rho\, e$.
The rest of the solver operates without any EOS branching.

### Constructor validation

| Provided | Result |
|---|---|
| `gamma_ideal_gas` only | Ideal-gas EOS |
| `C_0`, `s`, `Gamma_0` (all three) | Mie-Gruneisen EOS |
| `gamma_ideal_gas` + any of `C_0`/`s`/`Gamma_0` | `ValueError` |
| Partial `C_0`/`s`/`Gamma_0` (1 or 2 of 3) | `ValueError` |
| Neither | `ValueError` |

## References

- H.S. Udaykumar, L. Tran, D.M. Belk, K.J. Vanden,
  "An Eulerian method for computation of multimaterial impact with ENO
  shock-capturing and sharp interfaces,"
  *Journal of Computational Physics*, vol. 186, pp. 136–177, 2003.
  [doi:10.1016/S0021-9991(03)00027-5](https://doi.org/10.1016/S0021-9991(03)00027-5).

## Example Figures

### Aluminum piston — stress and velocity profiles

![Stress profile](assets/stress.png)

![Velocity profile](assets/velocity.png)

The three piecewise-constant regions are clearly visible: the **shocked
region** (behind the plastic shock), the **elastic region** (between the
plastic shock and the elastic precursor), and the **initial undisturbed
region** (ahead of the elastic precursor).

### High-strength case ($Y_0/G = 0.40$) — stress and velocity comparison

![High-strength comparison](assets/high_strength_comparison.png)

### Mie-Gruneisen EOS — all three regimes

Fictitious material: $\rho_0 = 4$, $C_0 = 7.5 \times 10^4$, $s = 1.3$,
$\Gamma_0 = 2$. Each plot overlays the original and energy-split
stress profiles.

**Elastic wave** ($G = 2 \times 10^{11}$, $Y_0 = 1.6 \times 10^{11}$,
$v_{piston} = 80000$):

![MG Elastic](assets/mg_elastic.png)

**Plastic shock and elastic wave** ($G = 2 \times 10^{11}$,
$Y_0 = 1 \times 10^{11}$, $v_{piston} = 90000$):

![MG Two-wave](assets/mg_two_wave.png)

**Plastic wave** ($G = 2 \times 10^{11}$, $Y_0 = 1 \times 10^{11}$,
$v_{piston} = 500000$):

![MG Plastic wave](assets/mg_plastic_wave.png)

### Ideal Gas EOS — all three regimes

Fictitious material: $\rho_0 = 2$, $\gamma = 4$, $G = 10^8$,
$Y_0 = 3 \times 10^7$, $e_{initial} = 0$.

**Elastic wave** ($v_{piston} = 1100$):

![IG Elastic](assets/ig_elastic.png)

**Plastic shock and elastic wave** ($v_{piston} = 2000$):

![IG Two-wave](assets/ig_two_wave.png)

**Plastic wave** ($v_{piston} = 10000$):

![IG Plastic wave](assets/ig_plastic_wave.png)

To regenerate all regime figures, run:

```bash
python examples/example_all_regimes.py
```