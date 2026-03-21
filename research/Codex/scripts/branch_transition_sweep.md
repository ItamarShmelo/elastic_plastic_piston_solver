# Strong-Shock Regime Investigation

## Setup

- Material family: `rho0=2.79`, `C0=5.33e5`, `s=1.34`, `Gamma0=2.0`, `G=2.86e11` (CGS).
- Map range: `Y0/G` from `0.01` to `1.0`, `v_piston` from `1e4` to `3e5 cm/s`.
- Classification rule requested: if any shock root exists above `U_se`, mark as `one_strong_shock`.

## Phase Map

![Regime phase map](../assets/branch_phase_map.png)

Grid resolution: `61 x 61` (`Y0/G` x `v_piston`), total `3721` cells per mode.

Status definitions:
- `three_wave`: roots only below `U_se`.
- `one_strong_shock`: at least one root above `U_se`.
- `no_plastic_shock`: `v_piston <= v_Y`.
- `no_root`: no detected sign-change root in scan range.

Coverage counts:
- Original: three_wave=129, one_strong_shock=1319, no_plastic_shock=2273, no_root=0
- Split: three_wave=119, one_strong_shock=1370, no_plastic_shock=2232, no_root=0

## Split-Mode Windows at Target Slices

At `v_piston = 80000 cm/s`:
- `three_wave` in `Y0/G`: [0.0100, 0.1276]
- `one_strong_shock` in `Y0/G`: [0.1338, 0.2266]

At `Y0/G = 0.15`:
- `three_wave` in `v_piston`: [51083, 77667] cm/s
- `one_strong_shock` in `v_piston`: [80083, 300000] cm/s

## Target Case Root Topology (`Y0/G=0.15`, `v_piston=80000`)

Original mode: `U_se=708600.94`, roots_below=[112580.34, 702771.23], roots_above=[], status=`three_wave`.
Split mode: `U_se=700411.66`, roots_below=[113287.05], roots_above=[701613.99], status=`one_strong_shock`.

## Reproduce

```bash
python3 research/Codex/scripts/example_branch_sweep.py
```

- Plot: `research/Codex/assets/branch_phase_map.png`
- Report: `research/Codex/scripts/branch_transition_sweep.md`
