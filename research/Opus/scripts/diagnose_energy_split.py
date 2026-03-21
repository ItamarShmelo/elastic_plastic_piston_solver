"""Diagnose the f_shock root landscape for the energy-split mode.

Evaluates f_shock(U_s) across the full range [v_piston+eps, U_se-eps]
for both the original and energy-split modes, finds all roots, and
reports physical quantities and constraint checks at each root.

Run:
    python examples/diagnose_energy_split.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

_PROJECT_ROOT = str(Path(__file__).resolve().parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

import numpy as np
from scipy.optimize import brentq

from elastoplastic_piston_solver import (
    _hugoniot_energy,
    _hugoniot_pressure,
    _mie_gruneisen_pressure,
)


def eos(e, rho, rho_0, C_0, s, Gamma_0):
    return _mie_gruneisen_pressure(e, rho, rho_0, C_0, s, Gamma_0)


def compute_elastic_state(params, energy_split):
    """Compute the elastic precursor state (rho_Y, P_Y, e_Y, U_se, v_Y)."""
    rho_0 = params["rho_0"]
    C_0 = params["C_0"]
    s = params["s"]
    Gamma_0 = params["Gamma_0"]
    G = params["G"]
    Y_0 = params["Y_0"]
    e_0 = params["e_initial"]

    two_thirds_Y0 = 2.0 / 3.0 * Y_0
    rho_Y = rho_0 * math.exp(Y_0 / (2.0 * G))
    e_el_Y = Y_0**2 / (6.0 * rho_Y * G)

    e_lo = e_0
    e_hi = e_0 + 10.0 * C_0**2

    if energy_split:
        def f_yield(e_th):
            P = eos(e_th, rho_Y, rho_0, C_0, s, Gamma_0)
            return (
                e_th + e_el_Y - e_0
                - 0.5 / (rho_Y * rho_0) * (P + two_thirds_Y0) * (rho_Y - rho_0)
            )
        e_th_Y = brentq(f_yield, e_lo, e_hi)
        e_Y = e_th_Y + e_el_Y
        P_Y = eos(e_th_Y, rho_Y, rho_0, C_0, s, Gamma_0)
    else:
        def f_yield(e):
            P = eos(e, rho_Y, rho_0, C_0, s, Gamma_0)
            return (
                e - e_0
                - 0.5 / (rho_Y * rho_0) * (P + two_thirds_Y0) * (rho_Y - rho_0)
            )
        e_Y = brentq(f_yield, e_lo, e_hi)
        e_th_Y = None
        e_el_Y_val = None
        P_Y = eos(e_Y, rho_Y, rho_0, C_0, s, Gamma_0)

    U_se_sq = rho_Y * (P_Y + two_thirds_Y0) / (rho_0 * (rho_Y - rho_0))
    U_se = math.sqrt(U_se_sq)
    v_Y = (rho_Y - rho_0) / rho_Y * U_se

    return {
        "rho_Y": rho_Y,
        "P_Y": P_Y,
        "e_Y": e_Y,
        "e_th_Y": e_th_Y,
        "e_el_Y": e_el_Y,
        "U_se": U_se,
        "v_Y": v_Y,
    }


def make_f_shock(params, elastic, energy_split):
    """Return the f_shock(U_s) function and a helper to get state at U_s."""
    rho_0 = params["rho_0"]
    C_0 = params["C_0"]
    s = params["s"]
    Gamma_0 = params["Gamma_0"]
    G = params["G"]
    Y_0 = params["Y_0"]
    v_piston = params["v_piston"]

    rho_Y = elastic["rho_Y"]
    P_Y = elastic["P_Y"]
    v_Y = elastic["v_Y"]
    four_thirds_Y0 = 4.0 / 3.0 * Y_0

    if energy_split:
        e_th_Y = elastic["e_th_Y"]
        e_el_Y = elastic["e_el_Y"]

        def f_shock(U_s):
            P_2 = P_Y + rho_Y * (U_s - v_Y) * (v_piston - v_Y)
            rho_2 = rho_Y * (U_s - v_Y) / (U_s - v_piston)
            e_el_2 = Y_0**2 / (6.0 * rho_2 * G)
            e_th_2 = (
                e_th_Y + e_el_Y
                + 0.5 / (rho_Y * rho_2) * (P_Y + P_2 + four_thirds_Y0) * (rho_2 - rho_Y)
                - e_el_2
            )
            return P_2 - eos(e_th_2, rho_2, rho_0, C_0, s, Gamma_0)

        def get_state(U_s):
            P_2 = P_Y + rho_Y * (U_s - v_Y) * (v_piston - v_Y)
            rho_2 = rho_Y * (U_s - v_Y) / (U_s - v_piston)
            e_el_2 = Y_0**2 / (6.0 * rho_2 * G)
            e_th_2 = (
                e_th_Y + e_el_Y
                + 0.5 / (rho_Y * rho_2) * (P_Y + P_2 + four_thirds_Y0) * (rho_2 - rho_Y)
                - e_el_2
            )
            return {
                "U_s": U_s,
                "P_2": P_2,
                "rho_2": rho_2,
                "e_total_2": e_th_2 + e_el_2,
                "e_th_2": e_th_2,
                "e_el_2": e_el_2,
                "rho_2/rho_Y": rho_2 / rho_Y,
                "mu_2": 1.0 - rho_0 / rho_2,
            }
    else:
        e_Y = elastic["e_Y"]

        def f_shock(U_s):
            P_2 = P_Y + rho_Y * (U_s - v_Y) * (v_piston - v_Y)
            rho_2 = rho_Y * (U_s - v_Y) / (U_s - v_piston)
            e_2 = (
                e_Y
                + 0.5 / (rho_Y * rho_2) * (P_Y + P_2 + four_thirds_Y0) * (rho_2 - rho_Y)
            )
            return P_2 - eos(e_2, rho_2, rho_0, C_0, s, Gamma_0)

        def get_state(U_s):
            P_2 = P_Y + rho_Y * (U_s - v_Y) * (v_piston - v_Y)
            rho_2 = rho_Y * (U_s - v_Y) / (U_s - v_piston)
            e_2 = (
                e_Y
                + 0.5 / (rho_Y * rho_2) * (P_Y + P_2 + four_thirds_Y0) * (rho_2 - rho_Y)
            )
            return {
                "U_s": U_s,
                "P_2": P_2,
                "rho_2": rho_2,
                "e_total_2": e_2,
                "e_th_2": None,
                "e_el_2": None,
                "rho_2/rho_Y": rho_2 / rho_Y,
                "mu_2": 1.0 - rho_0 / rho_2,
            }

    return f_shock, get_state


def find_all_roots(f, U_s_lo, U_s_hi, n_samples=5000):
    """Scan [U_s_lo, U_s_hi] and find all roots via sign changes."""
    Us_grid = np.linspace(U_s_lo, U_s_hi, n_samples)
    f_vals = np.array([f(u) for u in Us_grid])

    roots = []
    for i in range(len(Us_grid) - 1):
        if not (np.isfinite(f_vals[i]) and np.isfinite(f_vals[i + 1])):
            continue
        if f_vals[i] * f_vals[i + 1] < 0:
            root = brentq(f, Us_grid[i], Us_grid[i + 1])
            roots.append(root)
    return roots, Us_grid, f_vals


def check_physical(state, params, elastic):
    """Check physical constraints and return a list of violations."""
    violations = []
    s = params["s"]
    rho_0 = params["rho_0"]
    rho_Y = elastic["rho_Y"]
    U_se = elastic["U_se"]

    if state["rho_2"] <= rho_Y:
        violations.append(f"rho_2={state['rho_2']:.6g} <= rho_Y={rho_Y:.6g}")
    if state["U_s"] >= U_se:
        violations.append(f"U_s={state['U_s']:.6g} >= U_se={U_se:.6g}")
    if state["e_total_2"] <= 0:
        violations.append(f"e_total_2={state['e_total_2']:.6g} <= 0")
    if state["e_th_2"] is not None and state["e_th_2"] <= 0:
        violations.append(f"e_th_2={state['e_th_2']:.6g} <= 0")
    mu_singularity = 1.0 / s
    if state["mu_2"] >= mu_singularity:
        violations.append(
            f"mu_2={state['mu_2']:.6g} >= 1/s={mu_singularity:.6g} (Hugoniot singularity)"
        )
    return violations


def main():
    Y0_over_G = 0.15

    rho_0 = 2.79
    C_0 = 5.33e5
    s = 1.34
    Gamma_0 = 2.0
    cv = 9.0e6
    G = 2.86e11
    Y_0 = Y0_over_G * G
    v_piston = 80000.0

    params = {
        "rho_0": rho_0, "C_0": C_0, "s": s, "Gamma_0": Gamma_0,
        "G": G, "Y_0": Y_0, "e_initial": 0.0, "v_piston": v_piston,
    }

    print(f"Y_0/G = {Y0_over_G}")
    print(f"Y_0 = {Y_0:.6g} dyn/cm^2")
    print(f"v_piston = {v_piston} cm/s")
    print()

    for mode_name, split in [("Original (no split)", False), ("Energy-split", True)]:
        print("=" * 90)
        print(f"  Mode: {mode_name}")
        print("=" * 90)

        elastic = compute_elastic_state(params, energy_split=split)
        print(f"  rho_Y  = {elastic['rho_Y']:.8g}")
        print(f"  P_Y    = {elastic['P_Y']:.8g}")
        print(f"  e_Y    = {elastic['e_Y']:.8g}")
        if split:
            print(f"  e_th_Y = {elastic['e_th_Y']:.8g}")
            print(f"  e_el_Y = {elastic['e_el_Y']:.8g}")
        print(f"  U_se   = {elastic['U_se']:.8g}")
        print(f"  v_Y    = {elastic['v_Y']:.8g}")
        print()

        if v_piston <= elastic["v_Y"]:
            print(f"  *** v_piston <= v_Y: no plastic shock. Skipping.")
            print()
            continue

        f_shock, get_state = make_f_shock(params, elastic, energy_split=split)

        U_se = elastic["U_se"]
        eps = 1e-6
        U_s_lo = v_piston * (1.0 + eps)
        U_s_hi = U_se * (1.0 - eps)

        roots, Us_grid, f_vals = find_all_roots(f_shock, U_s_lo, U_s_hi, n_samples=10000)

        print(f"  f_shock range: [{U_s_lo:.6g}, {U_s_hi:.6g}]")
        print(f"  Number of roots found: {len(roots)}")
        print()

        for i, root in enumerate(roots):
            state = get_state(root)
            violations = check_physical(state, params, elastic)
            tag = "PHYSICAL" if not violations else "UNPHYSICAL"
            print(f"  Root {i+1}: U_s = {root:.8g}   [{tag}]")
            print(f"    rho_2      = {state['rho_2']:.8g}")
            print(f"    P_2        = {state['P_2']:.8g}")
            print(f"    e_total_2  = {state['e_total_2']:.8g}")
            if state["e_th_2"] is not None:
                print(f"    e_th_2     = {state['e_th_2']:.8g}")
                print(f"    e_el_2     = {state['e_el_2']:.8g}")
            print(f"    rho_2/rho_Y = {state['rho_2/rho_Y']:.8g}")
            print(f"    mu_2       = {state['mu_2']:.8g}")
            if violations:
                for v in violations:
                    print(f"    *** VIOLATION: {v}")
            print()

        # Show what the current solver would pick (first sign change from U_se)
        print("  --- Current solver bracket search (from U_se downward) ---")
        U_s_probe = U_s_hi
        f_at_top = f_shock(U_s_probe)
        ratio = 0.9
        bracket_lo = None
        for _ in range(200):
            U_s_probe = max(U_s_probe * ratio, v_piston * 1.001)
            f_probe = f_shock(U_s_probe)
            if f_at_top * f_probe < 0:
                bracket_lo = U_s_probe
                break
            if U_s_probe <= v_piston * 1.001:
                break

        if bracket_lo is not None:
            root_solver = brentq(f_shock, bracket_lo, U_s_hi)
            state_solver = get_state(root_solver)
            violations_solver = check_physical(state_solver, params, elastic)
            tag = "PHYSICAL" if not violations_solver else "UNPHYSICAL"
            print(f"  Solver picks U_s = {root_solver:.8g}   [{tag}]")
            print(f"    rho_2     = {state_solver['rho_2']:.8g}")
            print(f"    P_2       = {state_solver['P_2']:.8g}")
            print(f"    e_total_2 = {state_solver['e_total_2']:.8g}")
            if state_solver["e_th_2"] is not None:
                print(f"    e_th_2    = {state_solver['e_th_2']:.8g}")
            if violations_solver:
                for v in violations_solver:
                    print(f"    *** VIOLATION: {v}")
        else:
            print("  Solver FAILED to bracket (no sign change found).")
        print()

        # Plot
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

            finite_mask = np.isfinite(f_vals)
            ax1.plot(Us_grid[finite_mask], f_vals[finite_mask], linewidth=1.5)
            ax1.axhline(0, color="k", linewidth=0.5)
            for r in roots:
                ax1.axvline(r, color="red", linestyle="--", alpha=0.7, label=f"root @ {r:.4g}")
            ax1.set_xlabel("U_s (cm/s)")
            ax1.set_ylabel("f_shock(U_s)")
            ax1.set_title(f"f_shock landscape — {mode_name} (Y0/G={Y0_over_G})")
            ax1.legend(fontsize=9)
            ax1.grid(True, alpha=0.3)

            rho2_vals = np.array([
                elastic["rho_Y"] * (u - elastic["v_Y"]) / (u - v_piston)
                for u in Us_grid
            ])
            ax2.plot(Us_grid, rho2_vals, linewidth=1.5)
            ax2.axhline(elastic["rho_Y"], color="gray", linestyle=":", label="rho_Y")
            for r in roots:
                st = get_state(r)
                ax2.axvline(r, color="red", linestyle="--", alpha=0.7)
                ax2.plot(r, st["rho_2"], "ro", markersize=8)
            ax2.set_xlabel("U_s (cm/s)")
            ax2.set_ylabel("rho_2 (g/cm^3)")
            ax2.set_title(f"Compression at each U_s — {mode_name}")
            ax2.set_ylim(elastic["rho_Y"] * 0.95, elastic["rho_Y"] * 3.0)
            ax2.legend(fontsize=9)
            ax2.grid(True, alpha=0.3)

            fig.tight_layout()
            assets = Path(__file__).resolve().parent.parent / "assets"
            assets.mkdir(exist_ok=True)
            suffix = "split" if split else "orig"
            fig.savefig(assets / f"diag_fshock_{suffix}.png", dpi=150)
            plt.close(fig)
            print(f"  Plot saved: {assets / f'diag_fshock_{suffix}.png'}")
        except Exception as exc:
            print(f"  (Plotting skipped: {exc})")
        print()


if __name__ == "__main__":
    main()
