# Elastoplastic Piston Investigation (Codex Summary)

## Scope and Initial Trigger

This investigation started from the case:

- `Y0/G = 0.15`
- `v_piston = 80000 cm/s`
- material family: `rho0=2.79`, `C0=5.33e5`, `s=1.34`, `Gamma0=2.0`, `G=2.86e11` (CGS)

The observed issue was that **energy-split mode returned an unphysical branch** for this case, suggesting multiple mathematical shock roots with only one physically admissible.

The goals were to:

1. re-check derivation-to-implementation consistency,
2. explain why branch ambiguity appears,
3. map where physical three-wave solutions exist,
4. compare original vs split mode behavior across regimes,
5. identify representative and extreme three-wave parameter sets.

---

## Core Physics and Numerical Assumptions Tracked

The following equations/assumptions were explicitly tracked through the investigation:

- Mie-Gruneisen EOS is used in both modes.
- Energy split uses:
  - `e = e_th + e_el`
  - `e_el = S:S / (4 rho G)` and at yield `e_el^Y = Y0^2 / (6 rho^Y G)`
- Yield state is solved by root finding in:
  - original mode: total energy `e_Y`
  - split mode: thermal energy `e_th_Y` (with elastic part added in conservation balance)
- Plastic shock speed `U_s` is found from a nonlinear residual `f_shock(U_s)=0`.
- Physical three-wave ordering requires `U_s < U_se` (plastic shock behind elastic precursor).
- Regime classification (final user-requested definition):
  - `three_wave`: roots only below `U_se`
  - `one_strong_shock`: any root above `U_se`
  - `no_plastic_shock`: `v_piston <= v_Y`
  - `no_root`: no sign-change root found in scan range

---

## Code/Algorithm Changes Made During Investigation

The solver and analysis scripts were progressively hardened to make branch behavior explicit:

- Added robust root scanning over intervals (`_scan_sign_change_roots`) instead of relying on one local bracket.
- Refactored shock calculations so state construction and residual evaluation are separable (`shock_state` + `f_shock` flow).
- Added admissibility checks for candidate roots (finite values, monotonic compression, split thermal energy non-negative).
- Added diagnostics payload for root topology (`roots_below`, `roots_above`, etc.).
- Added explicit branch-breakdown failure behavior for cases where a valid moderate-compression branch is absent below `U_se`.
- Updated classification logic in phase-map tooling to match user interpretation:
  - **if any root exists above `U_se`, classify as `one_strong_shock`**.

---

## Chronological Plot Record (All Generated Plots, In Order)

### 1) Baseline stress profile (single-case aluminum)

![stress](assets/stress.png)

Used as an early baseline visualization of wave structure.

### 2) Baseline velocity profile (single-case aluminum)

![velocity](assets/velocity.png)

Companion baseline to confirm precursor/shock positioning.

### 3) Aluminum mode comparison (`Y0/G ~ 0.0091`)

![aluminum](assets/aluminum_comparison.png)

Original vs split nearly coincide (small elastic-energy fraction).

### 4) Copper mode comparison (`Y0/G ~ 0.002`)

![copper](assets/copper_comparison.png)

Again confirms negligible split/original differences at low `Y0/G`.

### 5) High-strength comparison (`Y0/G = 0.70`)

![high_strength](assets/high_strength_comparison.png)

Shows large divergence between modes (multi-percent to double-digit shifts).

### 6) High-velocity comparison at `Y0/G = 0.15`

![high_vp](assets/high_vp_mode_profile_comparison.png)

Used to inspect behavior at very large piston speed after branch-ambiguity discussions.

### 7) Phase map: three-wave vs strong-shock (original and split)

![phase_map](assets/branch_phase_map.png)

Final classification map over:

- `Y0/G in [0.01, 1.0]`
- `v_piston in [1e4, 3e5] cm/s`

Key split-mode windows from this report:

- at `v_piston = 80000`: three-wave for `Y0/G ~ [0.0100, 0.1276]`
- at `Y0/G = 0.15`: three-wave for `v_piston ~ [51083, 77667] cm/s`
- target case `Y0/G=0.15, v_piston=80000`:
  - original: roots below `U_se` only (`three_wave`)
  - split: one root below and one above `U_se` (`one_strong_shock`)

### 8) First zoom map (`U_se/U_s`) from earlier two-wave focus

![two_wave_zoom](assets/two_wave_zoom_use_over_us.png)

Early zoom before user correction of target region.

### 9) First zoom map, original vs split comparison

![two_wave_compare](assets/two_wave_zoom_use_over_us_compare_modes.png)

Companion panel for the same early zoom scope.

### 10) Corrected zoom on the green three-wave area (`U_se/U_s`)

![three_wave_zoom](assets/three_wave_zoom_use_over_us.png)

Updated after user clarification: zoom bounds follow the three-wave region.

### 11) Corrected zoom, original vs split (`U_se/U_s`)

![three_wave_compare](assets/three_wave_zoom_use_over_us_compare_modes.png)

Direct mode-to-mode view in the corrected three-wave window.

### 12) Full ratio atlas in three-wave overlap (`original/split`)

![ratio_maps](assets/three_wave_ratio_maps_orig_over_split.png)

Includes:

- `rho_Y`, `sigma_Y`, `v_Y`
- `rho_2`, `sigma_2`, `v_2`
- `U_se`, `U_s`

for cells where both modes remain three-wave.

### 13) Selected realistic three-wave case profile comparison

![selected_case](assets/selected_three_wave_case_comparison.png)

Map-guided case with clear but moderate split/original differences and distinct precursor.

### 14) Fictitious high-contrast three-wave case profile comparison

![fictitious_case](assets/fictitious_three_wave_case_comparison.png)

Final validated fictitious case constrained to be three-wave in **both** modes.

---

## Major Corrections and Iterative Decisions

The investigation included several important corrections driven by intermediate findings:

1. Initial branch interpretation was revised to user-required rule (`any root above U_se => one_strong_shock`).
2. Initial zoom targeted the wrong regime and was redone for the green three-wave region.
3. Initial fictitious candidate was rejected because it was not truly three-wave in both modes; search was rerun with explicit topology filtering.
4. Admissibility logic was adjusted so negative-energy rejection is applied specifically to split thermal energy constraints, not indiscriminately to total energy in original mode.

---

## Final Parameter Sets Highlighted

### Target problematic case (original question)

```python
Y0_over_G = 0.15
rho_0 = 2.79
C_0 = 5.33e5
s = 1.34
Gamma_0 = 2.0
G = 2.86e11
Y_0 = Y0_over_G * G
e_initial = 0.0
v_piston = 80000.0
```

### Selected three-wave case (Case A)

```python
rho_0 = 2.79
C_0 = 5.33e5
s = 1.34
Gamma_0 = 2.0
G = 2.86e11
Y_0 = 0.15352941176470586 * G
e_initial = 0.0
v_piston = 52666.66666666667
```

### Final validated fictitious three-wave case (Case B)

```python
rho_0 = 1.6522415952717888
C_0 = 288552.5906447974
s = 1.785505452715448
Gamma_0 = 3.2661667744752583
G = 7.12454606641788e11
Y_0 = 2.2003253037787656e11   # Y0/G = 0.308837
e_initial = 0.0
v_piston = 144962.60803487623

# Plot setup used
t_plot = 1.0e-6  # s
x_min = 0.0
x_max = 1.2019889680225462  # cm
x_limits = (x_min, x_max)
```

For this fictitious case, relative split-vs-original changes were large while staying three-wave:

- `U_s = -6.9840%`
- `U_se = -11.3673%`
- `rho_2 = +2.4428%`
- `P_2 = -31.6277%`
- `e_2 = +1.8724%`
- `v_Y = -11.3673%`

---

## Conclusions

1. The unphysical outcome at `Y0/G=0.15, v_piston=80000` is a **branch-selection/topology issue**, not a single algebraic typo.
2. Energy-split mode can introduce/shift branch structure so that a strong-shock root appears above `U_se` while a lower root still exists.
3. Robust root scanning plus explicit admissibility and topology diagnostics are required for reliable physical interpretation.
4. For low `Y0/G` (typical metals), split and original modes are nearly identical.
5. Differences grow rapidly with increasing `Y0/G`, and can become very large while still maintaining a three-wave structure (as shown by the validated fictitious case).

---

## Reproducibility (Packaged Scripts)

All scripts below are packaged in `research/Codex/scripts/` and import the local solver copy at `research/Codex/elastoplastic_piston_solver.py`.

```bash
python3 research/Codex/scripts/example_aluminum_piston.py
python3 research/Codex/scripts/example_mode_comparison.py
python3 research/Codex/scripts/example_high_vp_mode_profile.py
python3 research/Codex/scripts/example_branch_sweep.py
python3 research/Codex/scripts/example_two_wave_zoom_map.py
python3 research/Codex/scripts/example_three_wave_ratio_maps.py
python3 research/Codex/scripts/example_selected_three_wave_profiles.py
python3 research/Codex/scripts/example_branch_regression.py
```

Generated markdown artifacts copied into the same packaged scripts folder:

- `research/Codex/scripts/energy_split_comparison.md`
- `research/Codex/scripts/branch_transition_sweep.md`
- `research/Codex/scripts/selected_three_wave_cases.md`
