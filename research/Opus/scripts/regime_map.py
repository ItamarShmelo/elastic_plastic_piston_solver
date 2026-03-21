"""Generate a regime map of v_piston vs Y0/G for the energy-split mode.

Classifies every (Y0/G, v_piston) point into one of:
  - "Physical 2-wave": elastic precursor + weak plastic shock (the correct solution)
  - "Strong shock only": only the strong-shock root survives (unphysical)
  - "No plastic shock": v_piston < v_Y, purely elastic response
  - "No root": f_shock has no sign change at all
  - "Overdriven": v_piston >= U_se

Produces:
  - A color-map PNG saved to assets/
  - A comprehensive markdown report saved to examples/

Run:
    python examples/regime_map.py
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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

from elastoplastic_piston_solver import _mie_gruneisen_pressure

ASSETS_DIR = Path(__file__).resolve().parent.parent / "assets"
REPORT_PATH = Path(__file__).resolve().parent / "regime_map_report.md"

# ── Material parameters (Aluminum-like, CGS) ──────────────────────────
RHO_0 = 2.79
C_0 = 5.33e5
S_PARAM = 1.34
GAMMA_0 = 2.0
G = 2.86e11

# ── Regime codes ──────────────────────────────────────────────────────
NO_PLASTIC = 0
PHYSICAL_2WAVE = 1
STRONG_ONLY = 2
NO_ROOT = 3
OVERDRIVEN = 4

REGIME_LABELS = {
    NO_PLASTIC: "No plastic shock",
    PHYSICAL_2WAVE: "Physical 2-wave",
    STRONG_ONLY: "Strong shock only",
    NO_ROOT: "No root",
    OVERDRIVEN: "Overdriven",
}


GAMMA_IG = GAMMA_0 + 1  # 3.0 — matches MG thermal-pressure coefficient


def _eos_mg(e: float, rho: float) -> float:
    return _mie_gruneisen_pressure(e, rho, RHO_0, C_0, S_PARAM, GAMMA_0)


def _eos_ig(e: float, rho: float) -> float:
    return (GAMMA_IG - 1.0) * rho * e


_eos = _eos_mg  # default EOS (kept for backward compat)


def _elastic_state_split(Y0_over_G: float, eos_fn=None):
    """Compute elastic-precursor state for the energy-split mode."""
    eos_fn = eos_fn or _eos
    Y_0 = Y0_over_G * G
    rho_Y = RHO_0 * math.exp(Y_0 / (2.0 * G))
    e_el_Y = Y_0**2 / (6.0 * rho_Y * G)
    two_thirds_Y0 = 2.0 / 3.0 * Y_0

    def f_yield(e_th):
        P = eos_fn(e_th, rho_Y)
        return (
            e_th + e_el_Y
            - 0.5 / (rho_Y * RHO_0) * (P + two_thirds_Y0) * (rho_Y - RHO_0)
        )

    e_th_Y = brentq(f_yield, 0.0, 10.0 * C_0**2)
    P_Y = eos_fn(e_th_Y, rho_Y)
    U_se = math.sqrt(
        rho_Y * (P_Y + two_thirds_Y0) / (RHO_0 * (rho_Y - RHO_0))
    )
    v_Y = (rho_Y - RHO_0) / rho_Y * U_se
    return rho_Y, P_Y, e_th_Y, e_el_Y, U_se, v_Y, Y_0


def _elastic_state_orig(Y0_over_G: float, eos_fn=None):
    """Compute elastic-precursor state for the original mode."""
    eos_fn = eos_fn or _eos
    Y_0 = Y0_over_G * G
    rho_Y = RHO_0 * math.exp(Y_0 / (2.0 * G))
    two_thirds_Y0 = 2.0 / 3.0 * Y_0

    def f_yield(e):
        P = eos_fn(e, rho_Y)
        return (
            e
            - 0.5 / (rho_Y * RHO_0) * (P + two_thirds_Y0) * (rho_Y - RHO_0)
        )

    e_Y = brentq(f_yield, 0.0, 10.0 * C_0**2)
    P_Y = eos_fn(e_Y, rho_Y)
    U_se = math.sqrt(
        rho_Y * (P_Y + two_thirds_Y0) / (RHO_0 * (rho_Y - RHO_0))
    )
    v_Y = (rho_Y - RHO_0) / rho_Y * U_se
    return rho_Y, P_Y, e_Y, U_se, v_Y, Y_0


def classify_point(Y0_over_G: float, v_piston: float, split: bool,
                   eos_fn=None) -> int:
    """Return the regime code for one (Y0/G, v_piston) point."""
    eos_fn = eos_fn or _eos
    try:
        if split:
            rho_Y, P_Y, e_ref, e_el_Y, U_se, v_Y, Y_0 = _elastic_state_split(
                Y0_over_G, eos_fn
            )
        else:
            rho_Y, P_Y, e_ref, U_se, v_Y, Y_0 = _elastic_state_orig(
                Y0_over_G, eos_fn
            )
            e_el_Y = Y_0**2 / (6.0 * rho_Y * G)
    except Exception:
        return NO_ROOT

    if v_piston <= v_Y:
        return NO_PLASTIC
    if v_piston >= U_se:
        return OVERDRIVEN

    four_thirds_Y0 = 4.0 / 3.0 * Y_0

    if split:
        def f_shock(U_s):
            P_2 = P_Y + rho_Y * (U_s - v_Y) * (v_piston - v_Y)
            rho_2 = rho_Y * (U_s - v_Y) / (U_s - v_piston)
            e_th_2 = (
                e_ref + e_el_Y
                + 0.5 / (rho_Y * rho_2)
                * (P_Y + P_2 + four_thirds_Y0) * (rho_2 - rho_Y)
                - Y_0**2 / (6.0 * rho_2 * G)
            )
            return P_2 - eos_fn(e_th_2, rho_2)
    else:
        def f_shock(U_s):
            P_2 = P_Y + rho_Y * (U_s - v_Y) * (v_piston - v_Y)
            rho_2 = rho_Y * (U_s - v_Y) / (U_s - v_piston)
            e_2 = (
                e_ref
                + 0.5 / (rho_Y * rho_2)
                * (P_Y + P_2 + four_thirds_Y0) * (rho_2 - rho_Y)
            )
            return P_2 - eos_fn(e_2, rho_2)

    U_lo = v_piston * (1.0 + 1e-6)
    U_hi = U_se * (1.0 - 1e-6)
    n_probe = 400
    probes = np.linspace(U_lo, U_hi, n_probe)
    try:
        fv = np.array([f_shock(u) for u in probes])
    except Exception:
        return NO_ROOT

    sign_changes = []
    for i in range(n_probe - 1):
        if np.isfinite(fv[i]) and np.isfinite(fv[i + 1]):
            if fv[i] * fv[i + 1] < 0:
                sign_changes.append(i)

    if not sign_changes:
        return NO_ROOT

    f_top = fv[-1] if np.isfinite(fv[-1]) else fv[np.isfinite(fv)][-1] if np.any(np.isfinite(fv)) else 0.0
    if f_top > 0 and len(sign_changes) >= 1:
        return PHYSICAL_2WAVE

    return STRONG_ONLY


def build_regime_grid(yg_vals, vp_vals, split: bool, eos_fn=None):
    """Classify every point on a 2-D grid."""
    ny, nv = len(yg_vals), len(vp_vals)
    grid = np.empty((ny, nv), dtype=int)
    for iy, yg in enumerate(yg_vals):
        for iv, vp in enumerate(vp_vals):
            grid[iy, iv] = classify_point(yg, vp, split, eos_fn)
    return grid


def find_2wave_upper_boundary(yg_fine, split: bool, eos_fn=None):
    """For each Y0/G, find the max v_piston that still yields PHYSICAL_2WAVE.

    Binary searches between v_Y (lower) and a generous upper bound.
    Returns NaN where no physical 2-wave regime exists.
    """
    eos_fn = eos_fn or _eos
    boundary_vp = []
    for yg in yg_fine:
        try:
            if split:
                _, _, _, _, U_se, v_Y, _ = _elastic_state_split(yg, eos_fn)
            else:
                _, _, _, U_se, v_Y, _ = _elastic_state_orig(yg, eos_fn)
        except Exception:
            boundary_vp.append(np.nan)
            continue

        lo = v_Y * 1.01
        hi = U_se * 0.99

        if classify_point(yg, lo, split, eos_fn) != PHYSICAL_2WAVE:
            boundary_vp.append(np.nan)
            continue

        if classify_point(yg, hi, split, eos_fn) == PHYSICAL_2WAVE:
            boundary_vp.append(hi)
            continue

        for _ in range(30):
            mid = (lo + hi) / 2.0
            if classify_point(yg, mid, split, eos_fn) == PHYSICAL_2WAVE:
                lo = mid
            else:
                hi = mid
        boundary_vp.append((lo + hi) / 2.0)
    return np.array(boundary_vp)


def solve_Us(Y0_over_G: float, v_piston: float, split: bool, eos_fn=None):
    """Return (U_s, U_se) for the physical weak-shock root, or (NaN, NaN)."""
    eos_fn = eos_fn or _eos
    try:
        if split:
            rho_Y, P_Y, e_ref, e_el_Y, U_se, v_Y, Y_0 = _elastic_state_split(
                Y0_over_G, eos_fn)
        else:
            rho_Y, P_Y, e_ref, U_se, v_Y, Y_0 = _elastic_state_orig(
                Y0_over_G, eos_fn)
            e_el_Y = Y_0**2 / (6.0 * rho_Y * G)
    except Exception:
        return np.nan, np.nan

    if v_piston <= v_Y or v_piston >= U_se:
        return np.nan, U_se

    four_thirds_Y0 = 4.0 / 3.0 * Y_0

    if split:
        def f_shock(U_s):
            P_2 = P_Y + rho_Y * (U_s - v_Y) * (v_piston - v_Y)
            rho_2 = rho_Y * (U_s - v_Y) / (U_s - v_piston)
            e_th_2 = (
                e_ref + e_el_Y
                + 0.5 / (rho_Y * rho_2)
                * (P_Y + P_2 + four_thirds_Y0) * (rho_2 - rho_Y)
                - Y_0**2 / (6.0 * rho_2 * G)
            )
            return P_2 - eos_fn(e_th_2, rho_2)
    else:
        def f_shock(U_s):
            P_2 = P_Y + rho_Y * (U_s - v_Y) * (v_piston - v_Y)
            rho_2 = rho_Y * (U_s - v_Y) / (U_s - v_piston)
            e_2 = (
                e_ref
                + 0.5 / (rho_Y * rho_2)
                * (P_Y + P_2 + four_thirds_Y0) * (rho_2 - rho_Y)
            )
            return P_2 - eos_fn(e_2, rho_2)

    U_lo = v_piston * (1.0 + 1e-6)
    U_hi = U_se * (1.0 - 1e-6)
    n_probe = 400
    probes = np.linspace(U_lo, U_hi, n_probe)
    try:
        fv = np.array([f_shock(u) for u in probes])
    except Exception:
        return np.nan, U_se

    roots = []
    for i in range(n_probe - 1):
        if np.isfinite(fv[i]) and np.isfinite(fv[i + 1]) and fv[i] * fv[i + 1] < 0:
            try:
                roots.append(brentq(f_shock, probes[i], probes[i + 1]))
            except Exception:
                pass

    if not roots:
        return np.nan, U_se

    return max(roots), U_se


def build_ratio_grid(yg_vals, vp_vals, split: bool, eos_fn=None):
    """Build a 2-D grid of log10(U_se/U_s - 1) (NaN outside physical 2-wave)."""
    ny, nv = len(yg_vals), len(vp_vals)
    ratio = np.full((ny, nv), np.nan)
    for iy, yg in enumerate(yg_vals):
        for iv, vp in enumerate(vp_vals):
            if classify_point(yg, vp, split, eos_fn) != PHYSICAL_2WAVE:
                continue
            U_s, U_se = solve_Us(yg, vp, split, eos_fn)
            if np.isfinite(U_s) and U_s > 0:
                r = U_se / U_s - 1.0
                if r > 0:
                    ratio[iy, iv] = math.log10(r)
    return ratio


def plot_ratio_map(ratio, yg_vals, vp_vals, title, filename,
                   boundary_data=None, vmin=None, vmax=None):
    """Color-map of log10(U_se/U_s - 1) inside the physical 2-wave wedge."""
    fig, ax = plt.subplots(figsize=(11, 7.5))
    vp_km = np.array(vp_vals) / 1e3

    masked = np.ma.masked_invalid(ratio)
    if vmin is None:
        vmin = np.nanmin(ratio) if np.any(np.isfinite(ratio)) else -4.0
    if vmax is None:
        vmax = np.nanmax(ratio) if np.any(np.isfinite(ratio)) else 0.0

    im = ax.pcolormesh(
        vp_km, yg_vals, masked, cmap="viridis",
        vmin=vmin, vmax=vmax, shading="nearest",
    )
    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label(r"$\log_{10}(U_{se}/U_s - 1)$", fontsize=14)

    ax.set_facecolor("#d9d9d9")

    if boundary_data:
        for label, yg_fine, bnd_vp, color, ls in boundary_data:
            mask = np.isfinite(bnd_vp)
            ax.plot(
                bnd_vp[mask] / 1e3, np.array(yg_fine)[mask],
                color=color, linewidth=2.2, linestyle=ls, label=label,
            )

    ax.set_xlabel(r"$v_{\mathrm{piston}}$  ($10^3$ cm/s)", fontsize=14)
    ax.set_ylabel(r"$Y_0 / G$", fontsize=14)
    ax.set_title(title, fontsize=15, pad=12)
    ax.tick_params(labelsize=11)

    if boundary_data:
        ax.legend(loc="upper left", fontsize=10, framealpha=0.92, edgecolor="#ccc")

    fig.tight_layout()
    ASSETS_DIR.mkdir(exist_ok=True)
    fig.savefig(ASSETS_DIR / filename, dpi=180)
    plt.close(fig)
    print(f"  Plot saved: {ASSETS_DIR / filename}")


def plot_regime_map(grid, yg_vals, vp_vals, title, filename, boundary_data=None):
    """Render a color-map of the regime grid."""
    cmap = ListedColormap([
        "#3b528b",   # NO_PLASTIC  - blue
        "#21918c",   # PHYSICAL    - teal
        "#fde725",   # STRONG      - yellow
        "#e15759",   # NO_ROOT     - red
        "#440154",   # OVERDRIVEN  - dark purple
    ])
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5]
    norm = BoundaryNorm(bounds, cmap.N)

    fig, ax = plt.subplots(figsize=(11, 7.5))

    vp_km = np.array(vp_vals) / 1e3
    im = ax.pcolormesh(
        vp_km, yg_vals, grid, cmap=cmap, norm=norm, shading="nearest"
    )

    if boundary_data:
        for label, yg_fine, bnd_vp, color, ls in boundary_data:
            mask = np.isfinite(bnd_vp)
            ax.plot(
                bnd_vp[mask] / 1e3, np.array(yg_fine)[mask],
                color=color, linewidth=2.2, linestyle=ls, label=label,
            )

    ax.set_xlabel(r"$v_{\mathrm{piston}}$  ($10^3$ cm/s)", fontsize=14)
    ax.set_ylabel(r"$Y_0 / G$", fontsize=14)
    ax.set_title(title, fontsize=15, pad=12)
    ax.tick_params(labelsize=11)

    legend_patches = [
        Patch(facecolor="#3b528b", label="No plastic shock ($v_p < v_Y$)"),
        Patch(facecolor="#21918c", label="Physical 2-wave (weak shock)"),
        Patch(facecolor="#fde725", edgecolor="#888", label="Strong shock only (unphysical)"),
        Patch(facecolor="#e15759", label="No root found"),
        Patch(facecolor="#440154", label=r"Overdriven ($v_p \geq U_{se}$)"),
    ]
    if boundary_data:
        from matplotlib.lines import Line2D
        for label, _, _, color, ls in boundary_data:
            legend_patches.append(
                Line2D([0], [0], color=color, linewidth=2.2, linestyle=ls, label=label)
            )
    ax.legend(
        handles=legend_patches, loc="upper left", fontsize=10,
        framealpha=0.92, edgecolor="#ccc",
    )

    fig.tight_layout()
    ASSETS_DIR.mkdir(exist_ok=True)
    fig.savefig(ASSETS_DIR / filename, dpi=180)
    plt.close(fig)
    print(f"  Plot saved: {ASSETS_DIR / filename}")


def _regime_stats(grid, yg_vals, vp_vals):
    """Return per-regime point counts and area fractions."""
    total = grid.size
    stats = {}
    for code, label in REGIME_LABELS.items():
        count = int(np.sum(grid == code))
        stats[label] = (count, count / total * 100)
    return stats


def build_report(
    grid_split, grid_orig, yg_vals, vp_vals,
    yg_fine, bnd_split, bnd_orig, bnd_no_plastic,
    *,
    grid_ig_split=None, grid_ig_orig=None,
    bnd_ig_split=None, bnd_ig_orig=None,
    bnd_no_plastic_ig=None,
):
    """Build the comprehensive markdown report."""
    has_ig = grid_ig_split is not None

    s = []

    # ── Title & intro ─────────────────────────────────────────────────
    s.append("# Regime Map Report: Mie-Gruneisen vs Ideal Gas EOS\n\n")
    s.append(
        "This report maps the solution regimes of the 1-D elastoplastic piston "
        "problem as a function of the yield-to-shear ratio $Y_0/G$ and the piston "
        "velocity $v_{\\text{piston}}$, for both the **original** (single total-energy) "
        "and **energy-split** modes, using two equations of state:\n\n"
        "1. **Mie-Gruneisen** (MG): $P = P_H(\\mu) + \\Gamma_0 \\rho (e - e_H(\\mu))$ "
        "with Hugoniot reference curve from the linear $U_s = C_0 + s\\,u_p$ fit.\n"
        "2. **Ideal gas** (IG): $P = (\\gamma - 1)\\,\\rho\\,e$ with "
        "$\\gamma = \\Gamma_0 + 1 = 3$.\n\n"
    )

    # ── Material parameters ───────────────────────────────────────────
    s.append("## Material parameters (CGS)\n\n")
    s.append(
        "| Parameter | Symbol | Value |\n"
        "|-----------|--------|-------|\n"
        f"| Reference density | $\\rho_0$ | {RHO_0} g/cm$^3$ |\n"
        f"| Bulk sound speed | $C_0$ | {C_0:.3e} cm/s |\n"
        f"| Hugoniot slope | $s$ | {S_PARAM} |\n"
        f"| Gruneisen parameter | $\\Gamma_0$ | {GAMMA_0} |\n"
        f"| Shear modulus | $G$ | {G:.3e} dyn/cm$^2$ |\n"
        f"| Yield stress | $Y_0$ | variable ($Y_0/G = 0.005$--$0.80$) |\n\n"
    )

    # ── Regime definitions ────────────────────────────────────────────
    s.append("## Solution regimes\n\n")
    s.append(
        "The plastic-shock speed $U_s$ is determined by finding roots of\n\n"
        "$$f_{\\text{shock}}(U_s) = P_2(U_s) - P_{\\text{eos}}\\bigl(e_2(U_s),\\, "
        "\\rho_2(U_s)\\bigr) = 0$$\n\n"
        "in the interval $(v_{\\text{piston}},\\, U_{se})$.  Depending on the "
        "number and character of the roots, the solution falls into one of five regimes:\n\n"
    )
    s.append(
        "| # | Regime | Color | Description |\n"
        "|---|--------|-------|-------------|\n"
        "| 1 | **No plastic shock** | Blue | $v_p < v_Y$: piston velocity below "
        "elastic particle velocity; purely elastic response. |\n"
        "| 2 | **Physical 2-wave** | Teal | Two roots; the weak-shock root "
        "(largest $U_s$, nearest $U_{se}$) gives moderate compression. "
        "Elastic precursor + plastic shock coexist. |\n"
        "| 3 | **Strong shock only** | Yellow | Only the strong-shock root "
        "survives ($U_s / U_{se} \\sim 0.15$--$0.3$, $\\rho_2/\\rho^Y "
        "\\sim 1.9$--$2.0$). Unphysical for the piston problem. |\n"
        "| 4 | **No root** | Red | $f_{\\text{shock}}$ has no sign change; "
        "no self-consistent shock exists. |\n"
        "| 5 | **Overdriven** | Purple | $v_p \\geq U_{se}$: piston outruns "
        "the elastic precursor. |\n\n"
    )

    # ── Regime maps ───────────────────────────────────────────────────
    s.append("## Regime maps — Mie-Gruneisen EOS\n\n")
    s.append("### Energy-split mode\n\n")
    s.append("![Energy-split regime map](../assets/regime_map_split.png)\n\n")
    s.append("### Original mode\n\n")
    s.append("![Original regime map](../assets/regime_map_orig.png)\n\n")
    s.append("### Comparison (both boundaries overlaid)\n\n")
    s.append("![Combined regime map](../assets/regime_map_combined.png)\n\n")

    if has_ig:
        s.append("## Regime maps — Ideal Gas EOS ($\\gamma = 3$)\n\n")
        s.append("### Energy-split mode\n\n")
        s.append("![IG energy-split regime map](../assets/regime_map_ig_split.png)\n\n")
        s.append("### Original mode\n\n")
        s.append("![IG original regime map](../assets/regime_map_ig_orig.png)\n\n")
        s.append("### MG vs IG boundary comparison\n\n")
        s.append("![MG vs IG boundaries](../assets/regime_map_ig_vs_mg.png)\n\n")

    # ── Ratio maps ─────────────────────────────────────────────────────
    s.append("## $\\log_{10}(U_{se}/U_s - 1)$ inside the physical 2-wave wedge\n\n")
    s.append(
        "The quantity $U_{se}/U_s - 1$ measures the fractional speed excess of "
        "the elastic precursor over the plastic shock. Plotted on a $\\log_{10}$ "
        "scale: $-3$ means the two waves differ by 0.1%, $-1$ by 10%, "
        "and $0$ by 100% (factor of 2). Grey areas are outside the "
        "physical 2-wave regime.\n\n"
    )
    s.append("### 2x2 panel (all four cases)\n\n")
    s.append("![Ratio panel](../assets/ratio_panel.png)\n\n")
    s.append("### Individual maps\n\n")
    s.append("| | Energy-split | Original |\n")
    s.append("|---|---|---|\n")
    s.append(
        "| **MG** | ![](../assets/ratio_mg_split.png) "
        "| ![](../assets/ratio_mg_orig.png) |\n"
    )
    s.append(
        "| **IG** | ![](../assets/ratio_ig_split.png) "
        "| ![](../assets/ratio_ig_orig.png) |\n\n"
    )

    # ── Regime statistics ─────────────────────────────────────────────
    s.append("## Regime statistics\n\n")
    stats_split = _regime_stats(grid_split, yg_vals, vp_vals)
    stats_orig = _regime_stats(grid_orig, yg_vals, vp_vals)
    if has_ig:
        stats_ig_split = _regime_stats(grid_ig_split, yg_vals, vp_vals)
        stats_ig_orig = _regime_stats(grid_ig_orig, yg_vals, vp_vals)
        s.append(
            "| Regime | MG split (%) | MG orig (%) "
            "| IG split (%) | IG orig (%) |\n"
            "|--------|-------------|------------"
            "|-------------|-------------|\n"
        )
        for label in REGIME_LABELS.values():
            _, ps = stats_split[label]
            _, po = stats_orig[label]
            _, pis = stats_ig_split[label]
            _, pio = stats_ig_orig[label]
            s.append(f"| {label} | {ps:.1f} | {po:.1f} | {pis:.1f} | {pio:.1f} |\n")
    else:
        s.append(
            "| Regime | Energy-split (%) | Original (%) |\n"
            "|--------|------------------|--------------|\n"
        )
        for label in REGIME_LABELS.values():
            _, pct_s = stats_split[label]
            _, pct_o = stats_orig[label]
            s.append(f"| {label} | {pct_s:.1f} | {pct_o:.1f} |\n")
    s.append("\n")

    # ── Boundary analysis ─────────────────────────────────────────────
    s.append("## Critical boundary: Physical 2-wave vs Strong shock only\n\n")
    s.append(
        "The boundary between the physical 2-wave regime and the strong-shock-only "
        "regime is governed by the sign of $f_{\\text{shock}}$ near $U_{se}$.\n\n"
        "- When $f_{\\text{shock}}(U_{se}^-) > 0$, the function starts positive "
        "near $U_{se}$ and crosses zero twice (weak + strong root).\n"
        "- When $f_{\\text{shock}}(U_{se}^-) < 0$, the weak-shock root has "
        "vanished and only the strong-shock root remains.\n\n"
    )
    s.append(
        "At the small-$Y_0/G$ limit, the critical piston velocity is approximately\n\n"
        f"$$v_p^* \\approx {0.163 * C_0:.0f} \\;\\text{{cm/s}}"
        f"\\quad (v_p^*/C_0 \\approx 0.163)$$\n\n"
        "This threshold is an intrinsic property of the Mie-Gruneisen EOS "
        "($C_0$, $s$, $\\Gamma_0$), independent of strength.\n\n"
    )
    s.append(
        "The energy-split mode shifts this boundary to **lower** $v_{\\text{piston}}$ "
        "because the EOS sees only $e_{\\text{th}}$ instead of $e_{\\text{total}}$. "
        "This lowers $P^Y$, which reduces the Hugoniot momentum pressure $P_2$ "
        "near $U_{se}$ more than it reduces $P_{\\text{eos}}$, causing "
        "$f_{\\text{shock}}$ to flip sign.\n\n"
    )

    # ── Boundary table ────────────────────────────────────────────────
    s.append("### Boundary table: maximum $v_p$ for physical 2-wave solution\n\n")
    s.append(
        "| $Y_0/G$ | $v_p^*$ split (cm/s) | $v_p^*/C_0$ split "
        "| $v_p^*$ original (cm/s) | $v_p^*/C_0$ orig | Shift |\n"
        "|---------|---------------------|----------------"
        "|------------------------|-----------------|-------|\n"
    )
    step = max(1, len(yg_fine) // 20)
    for i, yg in enumerate(yg_fine):
        if i % step != 0 and i != len(yg_fine) - 1:
            continue
        vs = bnd_split[i]
        vo = bnd_orig[i]
        vs_str = f"{vs:.0f}" if np.isfinite(vs) else "--"
        vo_str = f"{vo:.0f}" if np.isfinite(vo) else "--"
        vsc = f"{vs / C_0:.4f}" if np.isfinite(vs) else "--"
        voc = f"{vo / C_0:.4f}" if np.isfinite(vo) else "--"
        if np.isfinite(vs) and np.isfinite(vo) and vo > 0:
            shift = f"{(vo - vs) / vo * 100:+.2f}%"
        else:
            shift = "--"
        s.append(
            f"| {yg:.3f} | {vs_str} | {vsc} "
            f"| {vo_str} | {voc} | {shift} |\n"
        )
    s.append("\n")

    # ── Key findings ──────────────────────────────────────────────────
    s.append("## Key findings — Mie-Gruneisen\n\n")
    s.append(
        "1. **The loss of the weak-shock root is tied to the Hugoniot reference "
        "curve** and its singularity at $\\mu = 1/s$. Both modes lose the "
        "physical solution above $v_p / C_0 \\approx 0.163$.\n\n"
    )
    s.append(
        "2. **Energy splitting shifts the critical boundary to lower "
        "$v_{\\text{piston}}$.**  At $Y_0/G = 0.15$ the shift is large enough "
        "that $v_p = 80{,}000$ cm/s falls outside the physical 2-wave region "
        "in energy-split mode but remains inside it in the original mode.\n\n"
    )
    s.append(
        "3. **The physical 2-wave window is a triangular wedge** bounded below "
        "by $v_Y(Y_0/G)$ (diagonal line) and above by $v_p^*(Y_0/G)$ (near-horizontal "
        "curve). As $Y_0/G$ increases, $v_Y$ rises and $v_p^*$ drops, squeezing "
        "the wedge closed.\n\n"
    )
    s.append(
        "4. **At extreme velocities ($v_p / C_0 \\gtrsim 0.65$), "
        "$f_{\\text{shock}}$ loses all roots**, creating the red \"no root\" region. "
        "This is caused by the MG Hugoniot singularity "
        "($\\mu \\to 1/s$).\n\n"
    )

    if has_ig:
        s.append("## Key findings — Ideal Gas EOS\n\n")
        s.append(
            "5. **The ideal gas EOS eliminates the Hugoniot singularity.** "
            "With $P = (\\gamma - 1)\\rho e$, there is no $(1 - s\\mu)^2$ "
            "denominator. The pressure is well-defined for all densities.\n\n"
        )
        s.append(
            "6. **The \"strong shock only\" (yellow) region vanishes entirely.** "
            "Without the MG cold-curve contribution ($P_H$), the "
            "$f_{\\text{shock}}$ function no longer develops the intermediate "
            "root structure that produces the unphysical strong-shock solution. "
            "The transition is directly from physical 2-wave to no-root.\n\n"
        )
        s.append(
            "7. **The \"no root\" region expands, not shrinks.** "
            "This is the key surprise. Without the cold pressure $P_H$, "
            "the EOS thermal pressure $(\\gamma-1)\\rho e$ exceeds the "
            "Rankine-Hugoniot momentum pressure $P_2$ across the entire "
            "$U_s$ range at moderate-to-high piston velocities. This is "
            "**not** a singularity artifact — the function is well-behaved "
            "(no NaN/inf) but simply never crosses zero. The two-wave "
            "structure genuinely cannot satisfy both the RH jump conditions "
            "and the ideal gas EOS simultaneously at those velocities.\n\n"
        )
        s.append(
            "8. **The physical 2-wave wedge is much larger** than with MG. "
            "The upper boundary ($v_p^*$) roughly doubles, extending to "
            "$v_p / C_0 \\approx 0.3$ instead of $\\approx 0.16$. "
            "This is because the ideal gas EOS lacks the stiff cold-curve "
            "pressure that flips $f_{\\text{shock}}$ negative near "
            "$U_{se}$ in the MG case.\n\n"
        )
        s.append(
            "9. **The \"Overdriven\" (purple) region becomes visible** in the "
            "bottom-right corner, where $v_p$ exceeds $U_{se}$ at high "
            "$Y_0/G$. With MG the elastic precursor speed $U_{se}$ is much "
            "higher (due to $P_H$), pushing this region beyond the plotted "
            "range.\n\n"
        )
        s.append(
            "10. **The energy-split effect persists** with the ideal gas EOS, "
            "because elastic energy subtraction is EOS-independent. "
            "The split-mode 2-wave wedge is slightly smaller than the "
            "original-mode wedge, consistent with the MG behavior.\n\n"
        )

    # ── Practical guidance ────────────────────────────────────────────
    s.append("## Practical guidance\n\n")
    s.append(
        "When using the energy-split formulation, ensure the piston velocity "
        "satisfies\n\n"
        "$$v_Y(Y_0/G) < v_{\\text{piston}} < v_p^*(Y_0/G)$$\n\n"
        "The $v_p^*$ boundary can be pre-computed for a given material. "
        "The solver now issues a `warnings.warn()` when only the strong-shock "
        "root is found, alerting users that the result is unphysical.\n\n"
    )
    if has_ig:
        s.append(
            "The ideal gas EOS roughly doubles the physical 2-wave velocity "
            "window by removing the cold-curve pressure that causes the "
            "MG weak-shock root to vanish early. However, it replaces the "
            "MG \"strong shock only\" failure with a \"no root\" failure at "
            "similar velocities. Neither EOS produces valid two-wave solutions "
            "above $v_p / C_0 \\approx 0.3$ (IG) or $\\approx 0.16$ (MG). "
            "The ideal gas comparison confirms that the strong-shock-only "
            "artifact is specific to the MG cold curve, while the loss of "
            "*any* solution at high $v_p$ is a structural property of the "
            "two-wave elastoplastic model.\n"
        )

    return "".join(s)


def main():
    ASSETS_DIR.mkdir(exist_ok=True)

    # ── Grid parameters ───────────────────────────────────────────────
    yg_vals = np.concatenate([
        np.arange(0.005, 0.05, 0.005),
        np.arange(0.05, 0.30, 0.01),
        np.arange(0.30, 0.81, 0.05),
    ])
    vp_vals = np.concatenate([
        np.arange(2000, 20000, 2000),
        np.arange(20000, 100000, 5000),
        np.arange(100000, 400001, 10000),
    ])

    print("Computing energy-split regime grid ...")
    grid_split = build_regime_grid(yg_vals, vp_vals, split=True)

    print("Computing original-mode regime grid ...")
    grid_orig = build_regime_grid(yg_vals, vp_vals, split=False)

    # ── Boundary curves ───────────────────────────────────────────────
    print("Computing boundary curves ...")
    yg_fine = np.linspace(0.005, 0.25, 50)

    bnd_split = find_2wave_upper_boundary(yg_fine, split=True)
    bnd_orig = find_2wave_upper_boundary(yg_fine, split=False)

    bnd_no_plastic = []
    for yg in yg_fine:
        try:
            _, _, _, _, _, v_Y, _ = _elastic_state_split(yg)
            bnd_no_plastic.append(v_Y)
        except Exception:
            bnd_no_plastic.append(np.nan)
    bnd_no_plastic = np.array(bnd_no_plastic)

    # ── Plots ─────────────────────────────────────────────────────────
    print("Generating plots ...")

    plot_regime_map(
        grid_split, yg_vals, vp_vals,
        title="Energy-Split Mode — Solution Regime Map",
        filename="regime_map_split.png",
        boundary_data=[
            (r"$v_p = v_Y$ boundary", yg_fine, bnd_no_plastic, "white", "--"),
            ("Physical / strong-shock boundary (split)", yg_fine, bnd_split, "#ff6361", "-"),
        ],
    )

    plot_regime_map(
        grid_orig, yg_vals, vp_vals,
        title="Original Mode — Solution Regime Map",
        filename="regime_map_orig.png",
        boundary_data=[
            (r"$v_p = v_Y$ boundary", yg_fine, bnd_no_plastic, "white", "--"),
            ("Physical / strong-shock boundary (orig)", yg_fine, bnd_orig, "#ff6361", "-"),
        ],
    )

    # Combined overlay: show split grid, overlay both boundaries
    plot_regime_map(
        grid_split, yg_vals, vp_vals,
        title="Energy-Split Regime Map with Original-Mode Boundary",
        filename="regime_map_combined.png",
        boundary_data=[
            (r"$v_p = v_Y$ boundary", yg_fine, bnd_no_plastic, "white", "--"),
            ("Boundary (energy-split)", yg_fine, bnd_split, "#ff6361", "-"),
            ("Boundary (original)", yg_fine, bnd_orig, "#58a4b0", "-."),
        ],
    )

    # ── Ideal-gas EOS grids ──────────────────────────────────────────
    print("Computing ideal-gas energy-split regime grid ...")
    grid_ig_split = build_regime_grid(yg_vals, vp_vals, split=True, eos_fn=_eos_ig)

    print("Computing ideal-gas original-mode regime grid ...")
    grid_ig_orig = build_regime_grid(yg_vals, vp_vals, split=False, eos_fn=_eos_ig)

    print("Computing ideal-gas boundary curves ...")
    bnd_ig_split = find_2wave_upper_boundary(yg_fine, split=True, eos_fn=_eos_ig)
    bnd_ig_orig = find_2wave_upper_boundary(yg_fine, split=False, eos_fn=_eos_ig)

    bnd_no_plastic_ig = []
    for yg in yg_fine:
        try:
            _, _, _, _, _, v_Y_ig, _ = _elastic_state_split(yg, _eos_ig)
            bnd_no_plastic_ig.append(v_Y_ig)
        except Exception:
            bnd_no_plastic_ig.append(np.nan)
    bnd_no_plastic_ig = np.array(bnd_no_plastic_ig)

    print("Generating ideal-gas plots ...")

    plot_regime_map(
        grid_ig_split, yg_vals, vp_vals,
        title=r"Ideal Gas ($\gamma=3$) Energy-Split — Regime Map",
        filename="regime_map_ig_split.png",
        boundary_data=[
            (r"$v_p = v_Y$ boundary", yg_fine, bnd_no_plastic_ig, "white", "--"),
            ("Physical / strong-shock boundary", yg_fine, bnd_ig_split, "#ff6361", "-"),
        ],
    )

    plot_regime_map(
        grid_ig_orig, yg_vals, vp_vals,
        title=r"Ideal Gas ($\gamma=3$) Original Mode — Regime Map",
        filename="regime_map_ig_orig.png",
        boundary_data=[
            (r"$v_p = v_Y$ boundary", yg_fine, bnd_no_plastic_ig, "white", "--"),
            ("Physical / strong-shock boundary", yg_fine, bnd_ig_orig, "#ff6361", "-"),
        ],
    )

    plot_regime_map(
        grid_ig_split, yg_vals, vp_vals,
        title=r"Ideal Gas ($\gamma=3$): MG vs IG Boundaries (energy-split)",
        filename="regime_map_ig_vs_mg.png",
        boundary_data=[
            (r"$v_p = v_Y$ (IG)", yg_fine, bnd_no_plastic_ig, "white", "--"),
            ("Boundary (IG, energy-split)", yg_fine, bnd_ig_split, "#ff6361", "-"),
            ("Boundary (MG, energy-split)", yg_fine, bnd_split, "#58a4b0", "-."),
        ],
    )

    # ── U_se / U_s ratio maps (zoomed to 2-wave wedge) ─────────────
    print("Computing U_se/U_s ratio grids (zoomed) ...")
    yg_zoom = np.linspace(0.005, 0.70, 70)
    vp_zoom = np.linspace(1000, 200000, 100)

    cases = [
        ("MG energy-split",  True,  _eos_mg, "ratio_mg_split.png",
         bnd_no_plastic, bnd_split),
        ("MG original",      False, _eos_mg, "ratio_mg_orig.png",
         bnd_no_plastic, bnd_orig),
        ("IG energy-split",  True,  _eos_ig, "ratio_ig_split.png",
         bnd_no_plastic_ig, bnd_ig_split),
        ("IG original",      False, _eos_ig, "ratio_ig_orig.png",
         bnd_no_plastic_ig, bnd_ig_orig),
    ]
    ratio_grids = {}
    for label, split, eos_fn, fname, bnd_vY, bnd_upper in cases:
        eos_tag = "IG" if eos_fn is _eos_ig else "MG"
        mode_tag = "energy-split" if split else "original"
        print(f"  {label} ...")
        ratio = build_ratio_grid(yg_zoom, vp_zoom, split, eos_fn)
        ratio_grids[label] = ratio

        bnd_vY_zoom = []
        for yg in yg_zoom:
            try:
                if split:
                    _, _, _, _, _, vY_z, _ = _elastic_state_split(yg, eos_fn)
                else:
                    _, _, _, _, vY_z, _ = _elastic_state_orig(yg, eos_fn)
                bnd_vY_zoom.append(vY_z)
            except Exception:
                bnd_vY_zoom.append(np.nan)
        bnd_vY_zoom = np.array(bnd_vY_zoom)

        bnd_upper_zoom = find_2wave_upper_boundary(yg_zoom, split, eos_fn)

        plot_ratio_map(
            ratio, yg_zoom, vp_zoom,
            title=f"$\\log_{{10}}(U_{{se}}/U_s - 1)$ — {eos_tag} {mode_tag}",
            filename=fname,
            boundary_data=[
                (r"$v_p = v_Y$", yg_zoom, bnd_vY_zoom, "white", "--"),
                (r"$v_p = v_p^*$", yg_zoom, bnd_upper_zoom, "white", "-"),
            ],
        )

    # ── 2×2 panel ────────────────────────────────────────────────────
    print("Generating 2×2 ratio panel ...")
    fig, axes = plt.subplots(2, 2, figsize=(16, 12), sharex=True, sharey=True)
    vp_km_zoom = np.array(vp_zoom) / 1e3
    panel_cases = [
        (axes[0, 0], "MG energy-split",  "MG energy-split"),
        (axes[0, 1], "MG original",      "MG original"),
        (axes[1, 0], "IG energy-split",  "IG energy-split"),
        (axes[1, 1], "IG original",      "IG original"),
    ]
    vmin_all = 0.0
    vmax_all = -10.0
    for _, key, _ in panel_cases:
        r = ratio_grids[key]
        if np.any(np.isfinite(r)):
            vmin_all = min(vmin_all, np.nanmin(r))
            vmax_all = max(vmax_all, np.nanmax(r))
    for ax, key, title_str in panel_cases:
        masked = np.ma.masked_invalid(ratio_grids[key])
        im = ax.pcolormesh(
            vp_km_zoom, yg_zoom, masked, cmap="viridis",
            vmin=vmin_all, vmax=vmax_all, shading="nearest",
        )
        ax.set_facecolor("#d9d9d9")
        ax.set_title(title_str, fontsize=13)
        ax.tick_params(labelsize=10)
    for ax in axes[1, :]:
        ax.set_xlabel(r"$v_{\mathrm{piston}}$  ($10^3$ cm/s)", fontsize=12)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$Y_0 / G$", fontsize=12)
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.025, 0.70])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(r"$\log_{10}(U_{se}/U_s - 1)$", fontsize=14)
    fig.suptitle(
        r"$\log_{10}(U_{se}/U_s - 1)$ in the Physical 2-Wave Regime",
        fontsize=15, y=0.97,
    )
    fig.savefig(ASSETS_DIR / "ratio_panel.png", dpi=180)
    plt.close(fig)
    print(f"  Plot saved: {ASSETS_DIR / 'ratio_panel.png'}")

    # ── Report ────────────────────────────────────────────────────────
    print("Writing report ...")
    report = build_report(
        grid_split, grid_orig, yg_vals, vp_vals,
        yg_fine, bnd_split, bnd_orig, bnd_no_plastic,
        grid_ig_split=grid_ig_split,
        grid_ig_orig=grid_ig_orig,
        bnd_ig_split=bnd_ig_split,
        bnd_ig_orig=bnd_ig_orig,
        bnd_no_plastic_ig=bnd_no_plastic_ig,
    )
    REPORT_PATH.write_text(report)
    print(f"  Report saved: {REPORT_PATH}")


if __name__ == "__main__":
    main()
