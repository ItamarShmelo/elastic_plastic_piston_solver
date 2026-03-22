"""Compare original and energy-split modes across three material cases.

All cases use CGS units.

Produces:
  - Overlay plots of stress and velocity profiles for a high-Y_0/G case
    where the two modes diverge visibly.
  - Quantitative comparison tables (printed) for:
      1. Aluminum  -- baseline, small Y_0/G
      2. Copper    -- baseline, small Y_0/G
      3. High-strength material -- Y_0/G ~ 0.70, large differences

Run:
    python examples/example_mode_comparison.py
"""

from __future__ import annotations

import sys
from pathlib import Path

_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

import re

import numpy as np
import matplotlib.pyplot as plt

from elastoplastic_piston_solver import ElastoplasticPistonSolver


def _sci_to_latex(s: str) -> str:
    """Convert Python scientific notation to LaTeX.

    '3.64661e+09' -> '3.64661 \\times 10^{9}'
    '1.3703e+07'  -> '1.3703 \\times 10^{7}'
    Leaves plain numbers (no 'e') unchanged.
    """
    return re.sub(
        r"([0-9.]+)e([+-])0*(\d+)",
        lambda m: (
            f"{m.group(1)} \\times 10^{{{'-' if m.group(2) == '-' else ''}{m.group(3)}}}"
        ),
        s,
    )

ASSETS_DIR = Path(__file__).resolve().parent.parent / "assets"

QUANTITIES = (
    "U_se", "U_s", "rho_Y", "P_Y", "v_Y", "e_Y",
    "rho_2", "P_2", "e_2",
)

QUANTITY_CONSOLE_LABELS = {
    "U_se": "U_se  (elastic precursor speed)",
    "U_s":  "U_s   (plastic shock speed)",
    "rho_Y": "rho^Y (yield density)",
    "P_Y":  "P^Y   (yield pressure)",
    "v_Y":  "v^Y   (yield particle velocity)",
    "e_Y":  "e^Y   (yield energy)",
    "rho_2": "rho_2 (shocked density)",
    "P_2":  "P_2   (shocked pressure)",
    "e_2":  "e_2   (shocked energy)",
}

QUANTITY_MD_LABELS = {
    "U_se":  "$U_{se}$",
    "U_s":   "$U_s$",
    "rho_Y": r"$\rho^Y$",
    "P_Y":   "$P^Y$",
    "v_Y":   "$v^Y$",
    "e_Y":   "$e^Y$",
    "rho_2": r"$\rho_2$",
    "P_2":   "$P_2$",
    "e_2":   "$e_2$",
}


def compare_solvers(
    label: str,
    kwargs: dict,
    unit_hint: str = "",
) -> tuple[ElastoplasticPistonSolver, ElastoplasticPistonSolver, list[tuple[str, float, float, float]]]:
    """Run both modes and print a comparison table."""
    s_orig = ElastoplasticPistonSolver(**kwargs, energy_split=False)
    s_split = ElastoplasticPistonSolver(**kwargs, energy_split=True)

    rows: list[tuple[str, float, float, float]] = []
    sep = "=" * 88
    print(sep)
    print(f"  {label}  {unit_hint}")
    print(sep)
    print(f"  {'Quantity':<40} {'Original':>16} {'Split':>16} {'Rel Diff':>10}")
    print("  " + "-" * 84)
    for name in QUANTITIES:
        v1 = getattr(s_orig, name)
        v2 = getattr(s_split, name)
        rd = abs(v2 - v1) / abs(v1) * 100 if abs(v1) > 0 else 0.0
        rows.append((name, v1, v2, rd))
        print(f"  {QUANTITY_CONSOLE_LABELS[name]:<40} {v1:>16.6g} {v2:>16.6g} {rd:>9.4f}%")

    if s_split.e_th_Y is not None and s_split.e_el_Y is not None:
        print()
        print(f"  {'Energy-split breakdown':<40} {'Thermal':>16} {'Elastic':>16}")
        print("  " + "-" * 84)
        print(f"  {'Yield region':<40} {s_split.e_th_Y:>16.6g} {s_split.e_el_Y:>16.6g}")
        assert s_split.e_th_2 is not None and s_split.e_el_2 is not None
        print(f"  {'Shocked region':<40} {s_split.e_th_2:>16.6g} {s_split.e_el_2:>16.6g}")

    print(sep)
    print()
    return s_orig, s_split, rows


def plot_comparison(
    s_orig: ElastoplasticPistonSolver,
    s_split: ElastoplasticPistonSolver,
    t: float,
    title_prefix: str,
    filename_prefix: str,
) -> None:
    """Produce stress and velocity overlay plots for both modes."""
    assert s_orig.U_se is not None and s_split.U_se is not None
    x_max = 1.2 * max(s_orig.U_se, s_split.U_se) * t
    x = np.linspace(0.0, x_max, 2000)

    r_orig = s_orig.solve(t, x)
    r_split = s_split.solve(t, x)

    LABEL_SIZE = 14
    TITLE_SIZE = 15
    LEGEND_SIZE = 12
    LW = 2.5

    fig, (ax_stress, ax_vel) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    PISTON_COLOR = "black"
    PISTON_LW = 1.8

    ax_stress.plot(x, r_orig["stress"], linewidth=LW, label="Original")
    ax_stress.plot(x, r_split["stress"], linewidth=LW, linestyle="--", label="Energy-split")
    ax_stress.axvline(r_orig["piston_location"], color=PISTON_COLOR,
                      linestyle="--", linewidth=PISTON_LW, label="Piston")
    ax_stress.set_ylabel(r"$\sigma_x = S_x - P$  (dyn/cm$^2$)", fontsize=LABEL_SIZE)
    ax_stress.set_title(f"{title_prefix} — Stress comparison", fontsize=TITLE_SIZE)
    ax_stress.legend(fontsize=LEGEND_SIZE)
    ax_stress.grid(True, alpha=0.3)
    ax_stress.tick_params(labelsize=11)

    ax_vel.plot(x, r_orig["velocity"], linewidth=LW, label="Original")
    ax_vel.plot(x, r_split["velocity"], linewidth=LW, linestyle="--", label="Energy-split")
    ax_vel.axvline(r_orig["piston_location"], color=PISTON_COLOR,
                   linestyle="--", linewidth=PISTON_LW)
    ax_vel.set_ylabel(r"$v$  (cm/s)", fontsize=LABEL_SIZE)
    ax_vel.set_xlabel(r"$x$  (cm)", fontsize=LABEL_SIZE)
    ax_vel.set_title(f"{title_prefix} — Velocity comparison", fontsize=TITLE_SIZE)
    ax_vel.legend(fontsize=LEGEND_SIZE)
    ax_vel.grid(True, alpha=0.3)
    ax_vel.tick_params(labelsize=11)

    fig.tight_layout()
    fig.savefig(ASSETS_DIR / f"{filename_prefix}_comparison.png", dpi=150)
    plt.close(fig)
    print(f"  Plot saved: {ASSETS_DIR / f'{filename_prefix}_comparison.png'}")


def build_markdown_table(
    rows: list[tuple[str, float, float, float]],
) -> str:
    """Return a markdown table string from comparison rows."""
    lines = [
        "| Quantity | Original | Energy-split | Rel. Diff (%) |",
        "|----------|----------|--------------|---------------|",
    ]
    for name, v1, v2, rd in rows:
        lines.append(
            f"| {QUANTITY_MD_LABELS[name]} "
            f"| ${_sci_to_latex(f'{v1:.6g}')}$ "
            f"| ${_sci_to_latex(f'{v2:.6g}')}$ "
            f"| ${rd:.4f}$ |"
        )
    return "\n".join(lines)


def main() -> None:
    ASSETS_DIR.mkdir(exist_ok=True)

    # ------------------------------------------------------------------
    # Case 1: Aluminum (CGS)
    # ------------------------------------------------------------------
    al_kw = dict(
        rho_0=2.79, C_0=5.33e5, s=1.34, Gamma_0=2.0,
        G=2.86e11, Y_0=2.6e9, e_initial=0.0, v_piston=5000.0,
    )
    s1_orig, s1_split, rows_al = compare_solvers(
        "Case 1: Aluminum", al_kw, "(CGS, Y_0/G = 0.0091)",
    )
    plot_comparison(s1_orig, s1_split, t=2.0e-6,
                    title_prefix="Aluminum", filename_prefix="aluminum")

    # ------------------------------------------------------------------
    # Case 2: Copper (CGS)
    # ------------------------------------------------------------------
    cu_kw = dict(
        rho_0=8.93, C_0=3.94e5, s=1.49, Gamma_0=2.0,
        G=4.5e11, Y_0=9.0e8, e_initial=0.0, v_piston=2000.0,
    )
    s2_orig, s2_split, rows_cu = compare_solvers(
        "Case 2: Copper", cu_kw, "(CGS, Y_0/G = 0.002)",
    )
    plot_comparison(s2_orig, s2_split, t=170.0e-6,
                    title_prefix="Copper", filename_prefix="copper")

    # ------------------------------------------------------------------
    # Case 3: High-strength material (CGS) — large Y_0/G
    # ------------------------------------------------------------------
    hs_kw = dict(
        rho_0=2.79, C_0=5.33e5, s=1.34, Gamma_0=2.0,
        G=2.86e11, Y_0=2.0e11, e_initial=0.0, v_piston=450000.0,
    )
    s3_orig, s3_split, rows_hs = compare_solvers(
        "Case 3: High-strength material", hs_kw, "(CGS, Y_0/G = 0.70)",
    )
    plot_comparison(s3_orig, s3_split, t=0.3e-6,
                    title_prefix=r"High-strength ($Y_0/G=0.70$)",
                    filename_prefix="high_strength")

    # ------------------------------------------------------------------
    # Write markdown report
    # ------------------------------------------------------------------
    md_path = Path(__file__).resolve().parent / "energy_split_comparison.md"
    md = _build_report(rows_al, rows_cu, rows_hs,
                       s1_split, s2_split, s3_split)
    md_path.write_text(md)
    print(f"  Report saved: {md_path}")


def _build_report(
    rows_al: list[tuple[str, float, float, float]],
    rows_cu: list[tuple[str, float, float, float]],
    rows_hs: list[tuple[str, float, float, float]],
    s_al: ElastoplasticPistonSolver,
    s_cu: ElastoplasticPistonSolver,
    s_hs: ElastoplasticPistonSolver,
) -> str:
    """Build the full markdown comparison report."""
    sections: list[str] = []

    sections.append("# Energy-Split Mode Comparison\n")
    sections.append(
        "This report compares the **original** (single total-energy) mode "
        "with the **energy-split** mode that decomposes internal energy into "
        "thermal and elastic parts: $e = e_{th} + e_{el}$, where "
        "$e_{el} = \\boldsymbol{S}:\\boldsymbol{S} \\,/\\, (4\\rho G)$.\n"
    )
    sections.append(
        "The EOS receives only the thermal energy $P(e_{th}, \\rho)$ "
        "in energy-split mode. The difference grows with the ratio $Y_0/G$.\n"
    )
    sections.append(
        "All cases use **CGS** units "
        "($\\mathrm{g/cm^3}$, $\\mathrm{cm/s}$, $\\mathrm{dyn/cm^2}$, $\\mathrm{erg/g}$).\n"
    )

    # Case 1 — Aluminum
    sections.append("## Case 1: Aluminum\n")
    sections.append(
        "Material parameters: "
        "$\\rho_0 = 2.79\\,\\mathrm{g/cm^3}$, "
        "$C_0 = 5.33 \\times 10^5\\,\\mathrm{cm/s}$, "
        "$s = 1.34$, "
        "$\\Gamma_0 = 2$, "
        "$G = 2.86 \\times 10^{11}\\,\\mathrm{dyn/cm^2}$, "
        "$Y_0 = 2.6 \\times 10^{9}\\,\\mathrm{dyn/cm^2}$, "
        "$v_{piston} = 5000\\,\\mathrm{cm/s}$.\n"
    )
    sections.append("$Y_0 / G = 0.0091$ — small elastic energy fraction.\n")
    sections.append(build_markdown_table(rows_al))
    sections.append("")
    assert s_al.e_th_Y is not None and s_al.e_el_Y is not None
    assert s_al.e_th_2 is not None and s_al.e_el_2 is not None
    sections.append(
        "Energy breakdown (split mode): "
        f"$e_{{th}}^Y = {_sci_to_latex(f'{s_al.e_th_Y:.6g}')}\\,\\mathrm{{erg/g}}$, "
        f"$e_{{el}}^Y = {_sci_to_latex(f'{s_al.e_el_Y:.6g}')}\\,\\mathrm{{erg/g}}$, "
        f"$e_{{th,2}} = {_sci_to_latex(f'{s_al.e_th_2:.6g}')}\\,\\mathrm{{erg/g}}$, "
        f"$e_{{el,2}} = {_sci_to_latex(f'{s_al.e_el_2:.6g}')}\\,\\mathrm{{erg/g}}$.\n"
    )
    sections.append("![Aluminum comparison](../assets/aluminum_comparison.png)\n")

    # Case 2 — Copper
    sections.append("## Case 2: Copper\n")
    sections.append(
        "Material parameters: "
        "$\\rho_0 = 8.93\\,\\mathrm{g/cm^3}$, "
        "$C_0 = 3.94 \\times 10^5\\,\\mathrm{cm/s}$, "
        "$s = 1.49$, "
        "$\\Gamma_0 = 2$, "
        "$G = 4.5 \\times 10^{11}\\,\\mathrm{dyn/cm^2}$, "
        "$Y_0 = 9.0 \\times 10^{8}\\,\\mathrm{dyn/cm^2}$, "
        "$v_{piston} = 2000\\,\\mathrm{cm/s}$.\n"
    )
    sections.append("$Y_0 / G = 0.002$ — very small elastic energy fraction.\n")
    sections.append(build_markdown_table(rows_cu))
    sections.append("")
    assert s_cu.e_th_Y is not None and s_cu.e_el_Y is not None
    assert s_cu.e_th_2 is not None and s_cu.e_el_2 is not None
    sections.append(
        "Energy breakdown (split mode): "
        f"$e_{{th}}^Y = {_sci_to_latex(f'{s_cu.e_th_Y:.6g}')}\\,\\mathrm{{erg/g}}$, "
        f"$e_{{el}}^Y = {_sci_to_latex(f'{s_cu.e_el_Y:.6g}')}\\,\\mathrm{{erg/g}}$, "
        f"$e_{{th,2}} = {_sci_to_latex(f'{s_cu.e_th_2:.6g}')}\\,\\mathrm{{erg/g}}$, "
        f"$e_{{el,2}} = {_sci_to_latex(f'{s_cu.e_el_2:.6g}')}\\,\\mathrm{{erg/g}}$.\n"
    )
    sections.append("![Copper comparison](../assets/copper_comparison.png)\n")

    # Case 3 — High-strength
    sections.append("## Case 3: High-Strength Material\n")
    sections.append(
        "To demonstrate a regime where the energy-split mode differs "
        "significantly, we use an artificial high-strength material with "
        "$Y_0 / G = 0.70$.\n"
    )
    sections.append(
        "Material parameters: "
        "$\\rho_0 = 2.79\\,\\mathrm{g/cm^3}$, "
        "$C_0 = 5.33 \\times 10^5\\,\\mathrm{cm/s}$, "
        "$s = 1.34$, "
        "$\\Gamma_0 = 2$, "
        "$G = 2.86 \\times 10^{11}\\,\\mathrm{dyn/cm^2}$, "
        "$Y_0 = 2.0 \\times 10^{11}\\,\\mathrm{dyn/cm^2}$, "
        "$v_{piston} = 4.5 \\times 10^5\\,\\mathrm{cm/s}$.\n"
    )
    sections.append(build_markdown_table(rows_hs))
    sections.append("")
    assert s_hs.e_th_Y is not None and s_hs.e_el_Y is not None
    assert s_hs.e_th_2 is not None and s_hs.e_el_2 is not None
    sections.append(
        "Energy breakdown (split mode): "
        f"$e_{{th}}^Y = {_sci_to_latex(f'{s_hs.e_th_Y:.6g}')}\\,\\mathrm{{erg/g}}$, "
        f"$e_{{el}}^Y = {_sci_to_latex(f'{s_hs.e_el_Y:.6g}')}\\,\\mathrm{{erg/g}}$, "
        f"$e_{{th,2}} = {_sci_to_latex(f'{s_hs.e_th_2:.6g}')}\\,\\mathrm{{erg/g}}$, "
        f"$e_{{el,2}} = {_sci_to_latex(f'{s_hs.e_el_2:.6g}')}\\,\\mathrm{{erg/g}}$.\n"
    )
    sections.append("![High-strength comparison](../assets/high_strength_comparison.png)\n")

    sections.append("## When Does the Energy-Split Mode Matter?\n")
    sections.append(
        "The elastic energy per unit mass at yield is "
        "$e_{el}^Y = Y_0^2 / (6 \\rho^Y G)$. "
        "Its ratio to the total internal energy determines how much "
        "the two modes diverge. For typical metals ($Y_0/G \\sim 10^{-3}$), "
        "the difference is negligible ($< 0.1\\%$). "
        "The effect becomes significant ($> 1\\%$) when $Y_0/G \\gtrsim 0.05$, "
        "which can occur in very high-strength or ceramic-like materials.\n"
    )

    return "\n".join(sections) + "\n"


if __name__ == "__main__":
    main()
