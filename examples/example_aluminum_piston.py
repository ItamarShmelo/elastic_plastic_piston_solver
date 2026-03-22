"""Example: Aluminum elastoplastic piston problem at t = 2 microseconds.

All values are in CGS units.

Material: Aluminum
    rho_0   = 2.79          g/cm^3
    C_0     = 5.33e5        cm/s
    s       = 1.34
    Gamma_0 = 2.0
    G       = 2.86e11       dyn/cm^2
    Y_0     = 2.6e9         dyn/cm^2

Initial conditions:
    e_initial = 0           erg/g
    v_piston  = 5000        cm/s
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


def main() -> None:
    # --- Material parameters (Aluminum, CGS) ---
    solver = ElastoplasticPistonSolver(
        rho_0=2.79,
        C_0=5.33e5,
        s=1.34,
        Gamma_0=2.0,
        G=2.86e11,
        Y_0=2.6e9,
        e_initial=0.0,
        v_piston=5000.0,
    )

    # --- Evaluate at t = 2 microseconds ---
    t: float = 2.0e-6  # seconds

    # Domain: from 0 to slightly beyond the elastic precursor
    assert solver.U_se is not None
    x_max: float = 1.2 * solver.U_se * t
    x = np.linspace(0.0, x_max, 1000)

    result = solver.solve(t, x)

    # --- Ensure output directory exists ---
    ASSETS_DIR.mkdir(exist_ok=True)

    # --- Plot styling constants ---
    LABEL_SIZE = 16
    TITLE_SIZE = 16
    LEGEND_SIZE = 13
    LINE_WIDTH = 2.5
    SHOCK_WIDTH = 1.8
    SHOCK_COLOR = "crimson"
    PRECURSOR_COLOR = "lime"
    t_label = r"$t = 2\;\mu\mathrm{s}$"

    # --- Combined plot ---
    fig, (ax_stress, ax_vel) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    ax_stress.plot(x, result["stress"], linewidth=LINE_WIDTH)
    ax_stress.set_ylabel(r"$\sigma_x = S_x - P$  (dyn/cm$^2$)", fontsize=LABEL_SIZE)
    ax_stress.set_title(r"$\sigma_x$ at " + t_label, fontsize=TITLE_SIZE)
    ax_stress.axvline(
        result["piston_location"],
        color="black",
        linestyle="--",
        linewidth=SHOCK_WIDTH,
        label="Piston",
    )
    ax_stress.axvline(
        result["shock_location"],
        color=SHOCK_COLOR,
        linestyle="--",
        linewidth=SHOCK_WIDTH,
        label="Plastic shock",
    )
    ax_stress.axvline(
        result["elastic_precursor_location"],
        color=PRECURSOR_COLOR,
        linestyle=":",
        linewidth=SHOCK_WIDTH,
        label="Elastic precursor",
    )
    ax_stress.legend(fontsize=LEGEND_SIZE)
    ax_stress.grid(True, alpha=0.3)
    ax_stress.tick_params(labelsize=12)

    ax_vel.plot(x, result["velocity"], linewidth=LINE_WIDTH)
    ax_vel.set_ylabel(r"$v$  (cm/s)", fontsize=LABEL_SIZE)
    ax_vel.set_xlabel(r"$x$  (cm)", fontsize=LABEL_SIZE)
    ax_vel.axvline(
        result["piston_location"],
        color="black",
        linestyle="--",
        linewidth=SHOCK_WIDTH,
    )
    ax_vel.axvline(
        result["shock_location"],
        color=SHOCK_COLOR,
        linestyle="--",
        linewidth=SHOCK_WIDTH,
    )
    ax_vel.axvline(
        result["elastic_precursor_location"],
        color=PRECURSOR_COLOR,
        linestyle=":",
        linewidth=SHOCK_WIDTH,
    )
    ax_vel.grid(True, alpha=0.3)
    ax_vel.tick_params(labelsize=12)

    fig.tight_layout()
    fig.savefig(ASSETS_DIR / "example_aluminum_piston.png", dpi=150)
    plt.close(fig)

    # --- Individual stress plot ---
    fig_stress, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x, result["stress"], linewidth=LINE_WIDTH)
    ax.set_ylabel(r"$\sigma_x = S_x - P$  (dyn/cm$^2$)", fontsize=LABEL_SIZE)
    ax.set_xlabel(r"$x$  (cm)", fontsize=LABEL_SIZE)
    ax.set_title(r"$\sigma_x$ at " + t_label, fontsize=TITLE_SIZE)
    ax.axvline(result["piston_location"], color="black", linestyle="--", linewidth=SHOCK_WIDTH, label="Piston")
    ax.axvline(result["shock_location"], color=SHOCK_COLOR, linestyle="--", linewidth=SHOCK_WIDTH, label="Plastic shock")
    ax.axvline(result["elastic_precursor_location"], color=PRECURSOR_COLOR, linestyle=":", linewidth=SHOCK_WIDTH, label="Elastic precursor")
    ax.legend(fontsize=LEGEND_SIZE)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=12)
    fig_stress.tight_layout()
    fig_stress.savefig(ASSETS_DIR / "stress.png", dpi=150)
    plt.close(fig_stress)

    # --- Individual velocity plot ---
    fig_vel, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x, result["velocity"], linewidth=LINE_WIDTH)
    ax.set_ylabel(r"$v$  (cm/s)", fontsize=LABEL_SIZE)
    ax.set_xlabel(r"$x$  (cm)", fontsize=LABEL_SIZE)
    ax.set_title(r"$v$ at " + t_label, fontsize=TITLE_SIZE)
    ax.axvline(result["piston_location"], color="black", linestyle="--", linewidth=SHOCK_WIDTH, label="Piston")
    ax.axvline(result["shock_location"], color=SHOCK_COLOR, linestyle="--", linewidth=SHOCK_WIDTH, label="Plastic shock")
    ax.axvline(result["elastic_precursor_location"], color=PRECURSOR_COLOR, linestyle=":", linewidth=SHOCK_WIDTH, label="Elastic precursor")
    ax.legend(fontsize=LEGEND_SIZE)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=12)
    fig_vel.tight_layout()
    fig_vel.savefig(ASSETS_DIR / "velocity.png", dpi=150)
    plt.close(fig_vel)

    print(f"Figures saved to {ASSETS_DIR}/")
    print(f"  Shock location:             {result['shock_location']:.6f} cm")
    print(f"  Elastic precursor location: {result['elastic_precursor_location']:.6f} cm")


if __name__ == "__main__":
    main()
