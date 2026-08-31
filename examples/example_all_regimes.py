"""Generate stress-profile figures for all 6 regime/EOS combinations.

Each figure overlays the original (no split) and energy-split stress
profiles with vertical dashed lines marking shock and precursor locations.

Mie-Gruneisen base: rho_0=4, C_0=7.5e4, s=1.3, Gamma_0=2.0
Ideal Gas:          rho_0=2, gamma=4, G=1e8, Y_0=3e7, e_initial=0

Run:
    python examples/example_all_regimes.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

import numpy as np
import matplotlib.pyplot as plt

from elastoplastic_piston_solver import ElastoplasticPistonSolver, WaveStructure

ASSETS_DIR = Path(__file__).resolve().parent.parent / "assets"

CASES: list[dict] = [
    dict(
        label="MG Elastic",
        filename="mg_elastic",
        params=dict(rho_0=4, C_0=7.5e4, s=1.3, Gamma_0=2.0,
                    G=2e11, Y_0=1.6e11, e_initial=0.0, v_piston=80000.0),
        expected=WaveStructure.ELASTIC_WAVE,
    ),
    dict(
        label="MG Two-wave",
        filename="mg_two_wave",
        params=dict(rho_0=4, C_0=7.5e4, s=1.3, Gamma_0=2.0,
                    G=2e11, Y_0=1e11, e_initial=0.0, v_piston=90000.0),
        expected=WaveStructure.PLASTIC_SHOCK_AND_ELASTIC_WAVE,
    ),
    dict(
        label="MG Plastic wave",
        filename="mg_plastic_wave",
        params=dict(rho_0=4, C_0=7.5e4, s=1.3, Gamma_0=2.0,
                    G=2e11, Y_0=1e11, e_initial=0.0, v_piston=500000.0),
        expected=WaveStructure.PLASTIC_WAVE,
    ),
    dict(
        label="IG Elastic",
        filename="ig_elastic",
        params=dict(rho_0=2, gamma_ideal_gas=4.0,
                    G=1e8, Y_0=3e7, e_initial=0.0, v_piston=1100.0),
        expected=WaveStructure.ELASTIC_WAVE,
    ),
    dict(
        label="IG Two-wave",
        filename="ig_two_wave",
        params=dict(rho_0=2, gamma_ideal_gas=4.0,
                    G=1e8, Y_0=3e7, e_initial=0.0, v_piston=2000.0),
        expected=WaveStructure.PLASTIC_SHOCK_AND_ELASTIC_WAVE,
    ),
    dict(
        label="IG Plastic wave",
        filename="ig_plastic_wave",
        params=dict(rho_0=2, gamma_ideal_gas=4.0,
                    G=1e8, Y_0=3e7, e_initial=0.0, v_piston=10000.0),
        expected=WaveStructure.PLASTIC_WAVE,
    ),
]

LABEL_SIZE = 14
TITLE_SIZE = 15
LEGEND_SIZE = 12
LW = 2.5
SHOCK_COLOR = "crimson"
PRECURSOR_COLOR = "lime"
SHOCK_STYLE = "--"
PRECURSOR_STYLE = ":"
MARKER_WIDTH = 1.8


def _reference_speed(s: ElastoplasticPistonSolver) -> float:
    """Return the fastest wave speed for domain sizing."""
    if s.wave_structure == WaveStructure.ELASTIC_WAVE:
        return s.U_se
    if s.wave_structure == WaveStructure.PLASTIC_WAVE:
        return s.U_s
    return s.U_se


def _plot_case(case: dict) -> None:
    params = case["params"]
    s_orig = ElastoplasticPistonSolver(**params, energy_split=False)
    s_split = ElastoplasticPistonSolver(**params, energy_split=True)

    assert s_orig.wave_structure == case["expected"], (
        f"{case['label']}: expected {case['expected'].name}, "
        f"got {s_orig.wave_structure.name}"
    )

    speed = max(_reference_speed(s_orig), _reference_speed(s_split))
    t = 1.0 / speed
    x = np.linspace(0.0, 1.3 * speed * t, 2000)

    r_orig = s_orig.solve(t, x)
    r_split = s_split.solve(t, x)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(x, r_orig["stress"], linewidth=LW, label="Original")
    ax.plot(x, r_split["stress"], linewidth=LW, linestyle="--",
            label="Energy-split")

    ax.axvline(r_orig["piston_location"], color="black", linestyle="--",
               linewidth=MARKER_WIDTH, label="Piston")

    for r, suffix in ((r_orig, ""), (r_split, " (split)")):
        loc_s = r["shock_location"]
        loc_p = r["elastic_precursor_location"]
        if math.isfinite(loc_s):
            ax.axvline(loc_s, color=SHOCK_COLOR, linestyle=SHOCK_STYLE,
                       linewidth=MARKER_WIDTH,
                       label=f"Plastic shock{suffix}" if not suffix else None)
        if math.isfinite(loc_p):
            ax.axvline(loc_p, color=PRECURSOR_COLOR, linestyle=PRECURSOR_STYLE,
                       linewidth=MARKER_WIDTH,
                       label=f"Elastic precursor{suffix}" if not suffix else None)

    vp = params["v_piston"]
    eos_tag = "Mie-Gruneisen" if "C_0" in params else "Ideal Gas"
    ax.set_title(
        f"{case['label']} — {eos_tag}, "
        rf"$v_{{piston}}={vp:.0f}$",
        fontsize=TITLE_SIZE,
    )
    ax.set_ylabel(r"$\sigma_x = S_x - P$", fontsize=LABEL_SIZE)
    ax.set_xlabel(r"$x$", fontsize=LABEL_SIZE)
    ax.legend(fontsize=LEGEND_SIZE)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=11)
    fig.tight_layout()
    fig.savefig(ASSETS_DIR / f"{case['filename']}.png", dpi=150)
    plt.close(fig)
    print(f"  {case['filename']}.png  ({s_orig.wave_structure.name})")


def main() -> None:
    ASSETS_DIR.mkdir(exist_ok=True)
    print("Generating regime figures:")
    for case in CASES:
        _plot_case(case)
    print(f"All figures saved to {ASSETS_DIR}/")


if __name__ == "__main__":
    main()
