"""Compare original vs energy-split modes using an Ideal Gas EOS.

Three cases are chosen to maximise the visible difference between modes:
  1. Y0/G = 0.20, v_piston = 40 000 cm/s  (60% P2 diff)
  2. Y0/G = 0.30, v_piston = 80 000 cm/s  (28% P2 diff, 5% U_s diff)
  3. Y0/G = 0.50, v_piston = 120 000 cm/s (36% P2 diff, 9% U_s diff)

Run:
    python examples/ig_mode_comparison.py
"""

from __future__ import annotations

import sys
from pathlib import Path

_PROJECT_ROOT = str(Path(__file__).resolve().parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

import warnings

import numpy as np
import matplotlib.pyplot as plt

from elastoplastic_piston_solver import ElastoplasticPistonSolver

ASSETS_DIR = Path(__file__).resolve().parent.parent / "assets"

# ── Material constants (CGS, Aluminum-like base) ──────────────────────
RHO_0 = 2.79
C_0 = 5.33e5
S_PARAM = 1.34
GAMMA_0 = 2.0
G = 2.86e11
GAMMA_IG = GAMMA_0 + 1  # 3.0


class IdealGasSolver(ElastoplasticPistonSolver):
    """ElastoplasticPistonSolver with an ideal-gas EOS P = (gamma-1)*rho*e."""

    def _eos(self, e: float, rho: float) -> float:
        return (GAMMA_IG - 1.0) * rho * e


CASES = [
    {
        "label": r"$Y_0/G = 0.20$, $v_p = 40\,000$ cm/s",
        "tag": "ig_case1",
        "Y0_over_G": 0.20,
        "v_piston": 40_000.0,
        "t": 1.5e-6,
    },
    {
        "label": r"$Y_0/G = 0.30$, $v_p = 80\,000$ cm/s",
        "tag": "ig_case2",
        "Y0_over_G": 0.30,
        "v_piston": 80_000.0,
        "t": 0.8e-6,
    },
    {
        "label": r"$Y_0/G = 0.50$, $v_p = 120\,000$ cm/s",
        "tag": "ig_case3",
        "Y0_over_G": 0.50,
        "v_piston": 120_000.0,
        "t": 0.5e-6,
    },
]


def make_solver(Y0_over_G: float, v_piston: float, split: bool):
    Y_0 = Y0_over_G * G
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return IdealGasSolver(
            rho_0=RHO_0, C_0=C_0, s=S_PARAM, Gamma_0=GAMMA_0,
            G=G, Y_0=Y_0, e_initial=0.0, v_piston=v_piston,
            energy_split=split,
        )


def run_case(case: dict):
    """Solve both modes, return (solver_orig, solver_split, result_orig, result_split)."""
    so = make_solver(case["Y0_over_G"], case["v_piston"], split=False)
    ss = make_solver(case["Y0_over_G"], case["v_piston"], split=True)

    t = case["t"]
    x_max = 1.25 * max(so.U_se, ss.U_se) * t
    x = np.linspace(0.0, x_max, 3000)

    ro = so.solve(t, x)
    rs = ss.solve(t, x)
    return so, ss, ro, rs, x


def pct(a, b):
    return abs(b - a) / abs(a) * 100 if abs(a) > 0 else 0.0


def plot_single_case(case, so, ss, ro, rs, x, filename):
    """3-panel plot: density, stress, velocity for one case."""
    LW = 2.2
    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True)

    # Density
    ax = axes[0]
    ax.plot(x, ro["density"], linewidth=LW, label="Original")
    ax.plot(x, rs["density"], linewidth=LW, linestyle="--", label="Energy-split")
    ax.set_ylabel(r"$\rho$  (g/cm$^3$)", fontsize=13)
    ax.set_title(f"Ideal Gas ($\\gamma={GAMMA_IG:.0f}$) — {case['label']}", fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # Stress
    ax = axes[1]
    ax.plot(x, ro["stress"], linewidth=LW, label="Original")
    ax.plot(x, rs["stress"], linewidth=LW, linestyle="--", label="Energy-split")
    ax.set_ylabel(r"$\sigma_x = S_x - P$  (dyn/cm$^2$)", fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # Velocity
    ax = axes[2]
    ax.plot(x, ro["velocity"], linewidth=LW, label="Original")
    ax.plot(x, rs["velocity"], linewidth=LW, linestyle="--", label="Energy-split")
    ax.set_ylabel(r"$v$  (cm/s)", fontsize=13)
    ax.set_xlabel(r"$x$  (cm)", fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    for a in axes:
        a.tick_params(labelsize=10)

    fig.tight_layout()
    ASSETS_DIR.mkdir(exist_ok=True)
    fig.savefig(ASSETS_DIR / filename, dpi=150)
    plt.close(fig)
    print(f"  Plot saved: {ASSETS_DIR / filename}")


def plot_combined_panel(all_data):
    """3×3 panel: rows = cases, cols = density / stress / velocity."""
    fig, axes = plt.subplots(3, 3, figsize=(18, 13))
    LW = 2.0

    col_labels = [
        r"$\rho$  (g/cm$^3$)",
        r"$\sigma_x$  (dyn/cm$^2$)",
        r"$v$  (cm/s)",
    ]
    field_keys = ["density", "stress", "velocity"]

    for irow, (case, so, ss, ro, rs, x) in enumerate(all_data):
        for icol, (key, ylabel) in enumerate(zip(field_keys, col_labels)):
            ax = axes[irow, icol]
            ax.plot(x, ro[key], linewidth=LW, label="Original")
            ax.plot(x, rs[key], linewidth=LW, linestyle="--", label="Energy-split")
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=9)

            if irow == 0:
                ax.set_title(ylabel, fontsize=13, pad=8)
            if icol == 0:
                ax.set_ylabel(case["label"], fontsize=11)
            if irow == len(all_data) - 1:
                ax.set_xlabel(r"$x$  (cm)", fontsize=11)
            if irow == 0 and icol == 0:
                ax.legend(fontsize=9, loc="best")

    fig.suptitle(
        f"Ideal Gas ($\\gamma={GAMMA_IG:.0f}$) — Original vs Energy-Split Profiles",
        fontsize=15, y=0.98,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    ASSETS_DIR.mkdir(exist_ok=True)
    fig.savefig(ASSETS_DIR / "ig_comparison_panel.png", dpi=150)
    plt.close(fig)
    print(f"  Plot saved: {ASSETS_DIR / 'ig_comparison_panel.png'}")


def print_table(case, so, ss):
    Y0_over_G = case["Y0_over_G"]
    sep = "=" * 90
    print(sep)
    print(f"  {case['label']}   (Y0/G = {Y0_over_G})")
    print(sep)
    print(f"  {'Quantity':<30} {'Original':>16} {'Split':>16} {'Rel Diff':>10}")
    print("  " + "-" * 76)
    for name, lbl in [
        ("U_se", "U_se (elastic speed)"),
        ("U_s",  "U_s  (plastic speed)"),
        ("v_Y",  "v_Y  (yield part. vel)"),
        ("rho_Y","rho_Y"),
        ("P_Y",  "P_Y"),
        ("rho_2","rho_2"),
        ("P_2",  "P_2"),
        ("e_2",  "e_2"),
    ]:
        vo = getattr(so, name)
        vs = getattr(ss, name)
        d = pct(vo, vs)
        print(f"  {lbl:<30} {vo:>16.6g} {vs:>16.6g} {d:>9.3f}%")
    print(f"  {'U_se / U_s':<30} {so.U_se/so.U_s:>16.6f} {ss.U_se/ss.U_s:>16.6f}")
    print(sep)
    print()


def main():
    ASSETS_DIR.mkdir(exist_ok=True)
    all_data = []

    for case in CASES:
        so, ss, ro, rs, x = run_case(case)
        print_table(case, so, ss)
        plot_single_case(case, so, ss, ro, rs, x, f"{case['tag']}_comparison.png")
        all_data.append((case, so, ss, ro, rs, x))

    plot_combined_panel(all_data)


if __name__ == "__main__":
    main()
