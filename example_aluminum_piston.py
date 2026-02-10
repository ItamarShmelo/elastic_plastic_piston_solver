"""Example: Aluminum elastoplastic piston problem at t = 2 microseconds.

All values are in CGS units.

Material: Aluminum
    rho_0   = 2.79          g/cm^3
    C_0     = 5.33e5        cm/s
    s       = 1.34
    Gamma_0 = 2.0
    P_0     = 0             dyn/cm^2
    G       = 2.86e11       dyn/cm^2
    Y_0     = 2.6e9         dyn/cm^2

Initial conditions:
    e_initial = 0           erg/g
    v_piston  = 5000        cm/s
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from elastoplastic_piston_solver import ElastoplasticPistonSolver


def main() -> None:
    # --- Material parameters (Aluminum, CGS) ---
    solver = ElastoplasticPistonSolver(
        rho_0=2.79,
        C_0=5.33e5,
        s=1.34,
        Gamma_0=2.0,
        P_0=0.0,
        G=2.86e11,
        Y_0=2.6e9,
        e_initial=0.0,
        v_piston=5000.0,
    )

    # --- Evaluate at t = 2 microseconds ---
    t: float = 2.0e-6  # seconds

    # Domain: from 0 to slightly beyond the elastic precursor
    x_max: float = 1.2 * solver.U_se * t
    x = np.linspace(0.0, x_max, 1000)

    result = solver.solve(t, x)

    # --- Plot ---
    fig, (ax_stress, ax_vel) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    ax_stress.plot(x, result["stress"], linewidth=1.5)
    ax_stress.set_ylabel(r"Total axial stress  $\sigma_x = S_x - P$  (dyn/cm$^2$)")
    ax_stress.set_title(f"Aluminum elastoplastic piston — t = {t * 1e6:.1f} µs")
    ax_stress.axvline(
        result["shock_location"],
        color="gray",
        linestyle="--",
        linewidth=0.8,
        label="Plastic shock",
    )
    ax_stress.axvline(
        result["elastic_precursor_location"],
        color="gray",
        linestyle=":",
        linewidth=0.8,
        label="Elastic precursor",
    )
    ax_stress.legend()
    ax_stress.grid(True, alpha=0.3)

    ax_vel.plot(x, result["velocity"], linewidth=1.5)
    ax_vel.set_ylabel("Velocity  $v$  (cm/s)")
    ax_vel.set_xlabel("Position  $x$  (cm)")
    ax_vel.axvline(
        result["shock_location"],
        color="gray",
        linestyle="--",
        linewidth=0.8,
    )
    ax_vel.axvline(
        result["elastic_precursor_location"],
        color="gray",
        linestyle=":",
        linewidth=0.8,
    )
    ax_vel.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig("example_aluminum_piston.png", dpi=150)
    print("Figure saved to example_aluminum_piston.png")
    print(f"  Shock location:             {result['shock_location']:.6f} cm")
    print(f"  Elastic precursor location: {result['elastic_precursor_location']:.6f} cm")


if __name__ == "__main__":
    main()
