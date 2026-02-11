"""Example: Copper elastoplastic piston problem — Validation against reference data.

All values are in SI units.

Material: Copper (OFHC)
    rho_0   = 8930          kg/m^3
    C_0     = 3940          m/s
    s       = 1.49
    Gamma_0 = 2.0
    G       = 45            GPa
    Y_0     = 90            MPa  (uniaxial yield stress)

Initial conditions:
    e_initial = 0           J/kg   (reference state for MG EOS)
    v_piston  = 20          m/s    (half the 40 m/s impact velocity)

The problem models symmetric impact: a body moving at 40 m/s strikes a
stationary body of the same material, producing a piston velocity of
20 m/s in the contact frame.

Reference
---------
Table 1: Exact values for variables for shocks of different strengths.
The elastoplastic case at 20 m/s piston velocity is validated:
  - Elastic precursor state
  - Plastic shock state
"""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure the project root is on the import path so the solver module is found
# regardless of where this script is invoked from.
_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

import numpy as np
import matplotlib.pyplot as plt

from elastoplastic_piston_solver import ElastoplasticPistonSolver

ASSETS_DIR = Path(__file__).resolve().parent.parent / "assets"


# ---------------------------------------------------------------------------
# Reference data from Table 1 (elastoplastic cases only)
# ---------------------------------------------------------------------------

REFERENCE = {
    "Elastic precursor, 20 m/s": {
        "wave_speed": 4722.0,
        "particle_velocity": 4.720,
        "density": 8938.9,
        "energy_jump": 11.138,
        "pressure": 139.03e6,
        "deviatoric_stress": -60.0e6,
        "total_stress": -199.03e6,
    },
    "Plastic shock, 20 m/s": {
        "wave_speed": 3977.0,
        "particle_velocity": 20.0,
        "density": 8973.5,
        "energy_jump": 213.53,
        "pressure": 681.59e6,
        "deviatoric_stress": -60.0e6,
        "total_stress": -741.59e6,
    },
}


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------

# Each entry: (display label, dict key, unit-conversion factor for display)
VARIABLES = [
    ("Wave speed (m/s)", "wave_speed", 1.0),
    ("Particle velocity (m/s)", "particle_velocity", 1.0),
    ("Density (kg/m³)", "density", 1.0),
    ("Internal energy jump (J/kg)", "energy_jump", 1.0),
    ("Pressure (MPa)", "pressure", 1.0e-6),
    ("Deviatoric stress (MPa)", "deviatoric_stress", 1.0e-6),
    ("Total stress (MPa)", "total_stress", 1.0e-6),
]


def print_comparison(
    computed: dict[str, dict[str, float]],
    reference: dict[str, dict[str, float]],
) -> None:
    """Print a formatted comparison table between computed and reference values."""
    sep = "=" * 90
    print(sep)
    print("Validation: Comparison against reference Table 1")
    print(sep)

    all_pass = True
    for case_name in computed:
        print(f"\n  {case_name}")
        print(f"  {'Variable':<35} {'Computed':>15} {'Reference':>15} {'Rel Error':>12}")
        print("  " + "-" * 78)
        for var_label, var_key, scale in VARIABLES:
            comp = computed[case_name][var_key] * scale
            ref = reference[case_name][var_key] * scale
            if abs(ref) > 0:
                rel_err = abs(comp - ref) / abs(ref) * 100
                err_str = f"{rel_err:.4f}%"
                if rel_err > 0.1:
                    all_pass = False
            else:
                err_str = "exact" if comp == 0.0 else f"|{abs(comp) * scale:.6f}|"
            print(f"  {var_label:<35} {comp:>15.4f} {ref:>15.4f} {err_str:>12}")

    print(f"\n{sep}")
    if all_pass:
        print("All computed values match reference data within 0.1% relative error.")
    else:
        print("WARNING: Some values exceed 0.1% relative error — check above.")
    print(sep)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    # --- Material parameters (Copper, SI) ---
    rho_0 = 8930.0       # kg/m^3
    C_0 = 3940.0         # m/s
    s = 1.49
    Gamma_0 = 2.0
    G = 45.0e9           # Pa
    Y_0 = 90.0e6         # Pa
    e_initial = 0.0      # J/kg  (MG EOS reference state)

    two_thirds_Y0 = 2.0 / 3.0 * Y_0  # deviatoric stress at yield

    # ------------------------------------------------------------------
    # 1. Elastoplastic solution for 20 m/s piston velocity
    # ------------------------------------------------------------------
    solver_20 = ElastoplasticPistonSolver(
        rho_0=rho_0,
        C_0=C_0,
        s=s,
        Gamma_0=Gamma_0,
        G=G,
        Y_0=Y_0,
        e_initial=e_initial,
        v_piston=20.0,
    )

    assert solver_20.U_se is not None and solver_20.v_Y is not None
    assert solver_20.rho_Y is not None and solver_20.e_Y is not None
    assert solver_20.P_Y is not None
    assert solver_20.U_s is not None and solver_20.rho_2 is not None
    assert solver_20.e_2 is not None and solver_20.P_2 is not None

    elastic_20: dict[str, float] = {
        "wave_speed": solver_20.U_se,
        "particle_velocity": solver_20.v_Y,
        "density": solver_20.rho_Y,
        "energy_jump": solver_20.e_Y - e_initial,
        "pressure": solver_20.P_Y,
        "deviatoric_stress": -two_thirds_Y0,
        "total_stress": -two_thirds_Y0 - solver_20.P_Y,
    }

    plastic_20: dict[str, float] = {
        "wave_speed": solver_20.U_s,
        "particle_velocity": solver_20.v_piston,
        "density": solver_20.rho_2,
        "energy_jump": solver_20.e_2 - e_initial,
        "pressure": solver_20.P_2,
        "deviatoric_stress": -two_thirds_Y0,
        "total_stress": -two_thirds_Y0 - solver_20.P_2,
    }

    # ------------------------------------------------------------------
    # 2. Compare computed values against reference Table 1
    # ------------------------------------------------------------------
    computed: dict[str, dict[str, float]] = {
        "Elastic precursor, 20 m/s": elastic_20,
        "Plastic shock, 20 m/s": plastic_20,
    }

    print_comparison(computed, REFERENCE)

    # ------------------------------------------------------------------
    # 3. Plot wave structure for the 20 m/s elastoplastic case
    # ------------------------------------------------------------------
    t: float = 170.0e-6  # 170 microseconds

    x_max: float = 1.2 * solver_20.U_se * t
    x = np.linspace(0.0, x_max, 2000)

    result = solver_20.solve(t, x)

    ASSETS_DIR.mkdir(exist_ok=True)

    # --- Plot styling ---
    LABEL_SIZE = 14
    TITLE_SIZE = 16
    LEGEND_SIZE = 11
    LINE_WIDTH = 2.5
    SHOCK_WIDTH = 1.8
    SHOCK_COLOR = "crimson"
    PRECURSOR_COLOR = "lime"
    t_label = r"$t = 170\;\mu\mathrm{s}$"

    fig, axes = plt.subplots(3, 2, figsize=(14, 12), sharex=True)

    def _decorate(ax: plt.Axes, ylabel: str, title: str = "") -> None:
        """Apply common styling to an axis."""
        ax.axvline(result["shock_location"], color=SHOCK_COLOR, ls="--",
                   lw=SHOCK_WIDTH, label="Plastic shock")
        ax.axvline(result["elastic_precursor_location"], color=PRECURSOR_COLOR,
                   ls=":", lw=SHOCK_WIDTH, label="Elastic precursor")
        ax.set_ylabel(ylabel, fontsize=LABEL_SIZE)
        if title:
            ax.set_title(title, fontsize=TITLE_SIZE)
        ax.legend(fontsize=LEGEND_SIZE)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=11)

    # (1,1) — Specific internal energy
    ax = axes[0, 0]
    ax.plot(x, result["energy"], linewidth=LINE_WIDTH)
    _decorate(ax, r"$e\;$ (J/kg)",
              title=f"Copper elastoplastic piston — {t_label}")

    # (1,2) — Pressure
    ax = axes[0, 1]
    ax.plot(x, np.array(result["pressure"]) * 1e-6, linewidth=LINE_WIDTH)
    _decorate(ax, r"$P\;$ (MPa)", title=f"{t_label}")

    # (2,1) — Deviatoric stress Sx
    ax = axes[1, 0]
    ax.plot(x, np.array(result["Sx"]) * 1e-6, linewidth=LINE_WIDTH)
    _decorate(ax, r"$S_x\;$ (MPa)")

    # (2,2) — Velocity
    ax = axes[1, 1]
    ax.plot(x, result["velocity"], linewidth=LINE_WIDTH)
    _decorate(ax, r"$v\;$ (m/s)")

    # (3,1) — Total stress sigma_x
    ax = axes[2, 0]
    ax.plot(x, np.array(result["stress"]) * 1e-6, linewidth=LINE_WIDTH)
    _decorate(ax, r"$\sigma_x = S_x - P\;$ (MPa)")
    ax.set_xlabel(r"$x\;$ (m)", fontsize=LABEL_SIZE)

    # (3,2) — Density
    ax = axes[2, 1]
    ax.plot(x, result["density"], linewidth=LINE_WIDTH)
    _decorate(ax, r"$\rho\;$ (kg/m$^3$)")
    ax.set_xlabel(r"$x\;$ (m)", fontsize=LABEL_SIZE)

    fig.tight_layout()
    fig.savefig(ASSETS_DIR / "example_copper_validation.png", dpi=150)
    plt.close(fig)

    print(f"\nFigure saved to {ASSETS_DIR / 'example_copper_validation.png'}")


if __name__ == "__main__":
    main()
