"""Generate high-v_piston profile comparison (original vs split).

Run:
    python3 research/Codex/scripts/example_high_vp_mode_profile.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from elastoplastic_piston_solver import ElastoplasticPistonSolver


def main() -> None:
    # High-v_piston case used in the investigation.
    Y0_over_G = 0.15
    rho0 = 2.79
    C0 = 5.33e5
    s = 1.34
    Gamma0 = 2.0
    G = 2.86e11
    Y0 = Y0_over_G * G
    v_piston = 350000.0

    common = dict(
        rho_0=rho0,
        C_0=C0,
        s=s,
        Gamma_0=Gamma0,
        G=G,
        Y_0=Y0,
        e_initial=0.0,
        v_piston=v_piston,
    )

    s_orig = ElastoplasticPistonSolver(**common, energy_split=False)
    s_split = ElastoplasticPistonSolver(**common, energy_split=True)

    t = 0.30e-6
    x_max = 1.20 * max(s_orig.U_se, s_split.U_se) * t
    x = np.linspace(0.0, x_max, 2400)

    r_orig = s_orig.solve(t, x)
    r_split = s_split.solve(t, x)

    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

    axes[0].plot(x, r_orig["stress"], lw=2.3, label="Original")
    axes[0].plot(x, r_split["stress"], lw=2.3, ls="--", label="Energy-split")
    axes[0].set_ylabel(r"$\sigma_x$ (dyn/cm$^2$)")
    axes[0].set_title(rf"High-$v_{{piston}}$ Comparison ($Y_0/G={Y0_over_G:.2f}$, $v_{{piston}}={v_piston:.0f}$ cm/s)")
    axes[0].grid(alpha=0.30)
    axes[0].legend()

    axes[1].plot(x, r_orig["velocity"], lw=2.3, label="Original")
    axes[1].plot(x, r_split["velocity"], lw=2.3, ls="--", label="Energy-split")
    axes[1].set_ylabel(r"$v$ (cm/s)")
    axes[1].grid(alpha=0.30)

    axes[2].plot(x, r_orig["density"], lw=2.3, label="Original")
    axes[2].plot(x, r_split["density"], lw=2.3, ls="--", label="Energy-split")
    axes[2].set_ylabel(r"$\rho$ (g/cm$^3$)")
    axes[2].set_xlabel("x (cm)")
    axes[2].grid(alpha=0.30)

    for ax in axes:
        ax.axvline(s_orig.U_s * t, color="tab:blue", alpha=0.25, lw=1.4)
        ax.axvline(s_orig.U_se * t, color="tab:blue", alpha=0.25, lw=1.4)
        ax.axvline(s_split.U_s * t, color="tab:orange", alpha=0.25, lw=1.4, ls="--")
        ax.axvline(s_split.U_se * t, color="tab:orange", alpha=0.25, lw=1.4, ls="--")

    fig.tight_layout()

    assets_dir = Path(__file__).resolve().parent.parent / "assets"
    assets_dir.mkdir(exist_ok=True)
    out_path = assets_dir / "high_vp_mode_profile_comparison.png"
    fig.savefig(out_path, dpi=180)
    plt.close(fig)

    print(f"Saved: {out_path}")
    print(f"Original: U_s={s_orig.U_s:.3f}, U_se={s_orig.U_se:.3f}, rho2={s_orig.rho_2:.5f}")
    print(f"Split:    U_s={s_split.U_s:.3f}, U_se={s_split.U_se:.3f}, rho2={s_split.rho_2:.5f}")


if __name__ == "__main__":
    main()
