"""Generate selected three-wave profile comparisons (original vs split).

Creates two cases:
1) A map-guided case with clear mode differences and a distinct elastic precursor.
2) A fictitious three-wave case with very large mode differences.

Run:
    python3 research/Codex/scripts/example_selected_three_wave_profiles.py
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
from example_branch_sweep import STATUS_THREE_WAVE, _classify_topology


def _rel(new: float, old: float) -> float:
    return (new - old) / old


def _plot_case(
    label: str,
    kwargs: dict[str, float],
    t: float,
    out_path: Path,
) -> dict[str, float]:
    s_orig = ElastoplasticPistonSolver(**kwargs, energy_split=False)
    s_split = ElastoplasticPistonSolver(**kwargs, energy_split=True)

    x_max = 1.2 * max(s_orig.U_se, s_split.U_se) * t
    x = np.linspace(0.0, x_max, 2600)
    r_orig = s_orig.solve(t, x)
    r_split = s_split.solve(t, x)

    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

    axes[0].plot(x, r_orig["density"], lw=2.3, label="Original")
    axes[0].plot(x, r_split["density"], lw=2.3, ls="--", label="Energy-split")
    axes[0].set_ylabel(r"$\rho$ (g/cm$^3$)")
    axes[0].set_title(label)
    axes[0].grid(alpha=0.3)
    axes[0].legend()

    axes[1].plot(x, r_orig["stress"], lw=2.3)
    axes[1].plot(x, r_split["stress"], lw=2.3, ls="--")
    axes[1].set_ylabel(r"$\sigma_x$ (dyn/cm$^2$)")
    axes[1].grid(alpha=0.3)

    axes[2].plot(x, r_orig["velocity"], lw=2.3)
    axes[2].plot(x, r_split["velocity"], lw=2.3, ls="--")
    axes[2].set_ylabel(r"$v$ (cm/s)")
    axes[2].set_xlabel("x (cm)")
    axes[2].grid(alpha=0.3)

    for ax in axes:
        ax.axvline(s_orig.U_s * t, color="tab:blue", alpha=0.28, lw=1.3)
        ax.axvline(s_orig.U_se * t, color="tab:blue", alpha=0.28, lw=1.3)
        ax.axvline(s_split.U_s * t, color="tab:orange", alpha=0.28, lw=1.3, ls="--")
        ax.axvline(s_split.U_se * t, color="tab:orange", alpha=0.28, lw=1.3, ls="--")

    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)

    return {
        "U_s_orig": s_orig.U_s,
        "U_s_split": s_split.U_s,
        "U_se_orig": s_orig.U_se,
        "U_se_split": s_split.U_se,
        "rho2_orig": s_orig.rho_2,
        "rho2_split": s_split.rho_2,
        "P2_orig": s_orig.P_2,
        "P2_split": s_split.P_2,
        "e2_orig": s_orig.e_2,
        "e2_split": s_split.e_2,
        "vY_orig": s_orig.v_Y,
        "vY_split": s_split.v_Y,
        "sep_orig": s_orig.U_se / s_orig.U_s,
        "sep_split": s_split.U_se / s_split.U_s,
        "t": t,
        "x_max": x_max,
    }


def _build_summary_md(path: Path, entries: list[dict[str, object]]) -> None:
    lines: list[str] = []
    lines.append("# Selected Three-Wave Case Comparisons")
    lines.append("")
    lines.append("All ratios below are `split/original - 1` (fraction).")
    lines.append("")

    for entry in entries:
        label = str(entry["label"])
        params = entry["params"]
        metrics = entry["metrics"]
        image = entry["image"]
        ratio = float(params["Y_0"]) / float(params["G"])

        lines.append(f"## {label}")
        lines.append("")
        lines.append(
            f"- Parameters: `Y0/G={ratio:.6f}`, `v_piston={params['v_piston']:.3f}` cm/s, "
            f"`rho0={params['rho_0']:.6g}`, `C0={params['C_0']:.6g}`, `s={params['s']:.6g}`, "
            f"`Gamma0={params['Gamma_0']:.6g}`, `G={params['G']:.6g}`"
        )
        lines.append(
            f"- Precursor separation: `U_se/U_s` original={metrics['sep_orig']:.6f}, "
            f"split={metrics['sep_split']:.6f}"
        )
        lines.append(
            f"- Relative deltas (split vs original): "
            f"`U_s={_rel(metrics['U_s_split'], metrics['U_s_orig']):+.4%}`, "
            f"`U_se={_rel(metrics['U_se_split'], metrics['U_se_orig']):+.4%}`, "
            f"`rho_2={_rel(metrics['rho2_split'], metrics['rho2_orig']):+.4%}`, "
            f"`P_2={_rel(metrics['P2_split'], metrics['P2_orig']):+.4%}`, "
            f"`e_2={_rel(metrics['e2_split'], metrics['e2_orig']):+.4%}`, "
            f"`v_Y={_rel(metrics['vY_split'], metrics['vY_orig']):+.4%}`"
        )
        lines.append(f"- Plot: ![]({image})")
        lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    assets = Path(__file__).resolve().parent.parent / "assets"
    assets.mkdir(exist_ok=True)

    # Case A (map-guided): clear mode difference + distinct elastic precursor.
    case_a = dict(
        rho_0=2.79,
        C_0=5.33e5,
        s=1.34,
        Gamma_0=2.0,
        G=2.86e11,
        Y_0=0.15352941176470586 * 2.86e11,
        e_initial=0.0,
        v_piston=52666.66666666667,
    )

    # Case B (fictitious): intentionally large split/original separation,
    # constrained to be three-wave in both modes.
    case_b = dict(
        rho_0=1.6522415952717888,
        C_0=288552.5906447974,
        s=1.785505452715448,
        Gamma_0=3.2661667744752583,
        G=712454606641.788,
        Y_0=220032530377.87656,
        e_initial=0.0,
        v_piston=144962.60803487623,
    )

    base_b = dict(
        rho_0=case_b["rho_0"],
        C_0=case_b["C_0"],
        s=case_b["s"],
        Gamma_0=case_b["Gamma_0"],
        G=case_b["G"],
        e_initial=case_b["e_initial"],
    )
    ratio_b = case_b["Y_0"] / case_b["G"]
    status_b_orig = int(
        _classify_topology(
            base=base_b,
            ratio=ratio_b,
            v_piston=case_b["v_piston"],
            energy_split=False,
            root_samples=4000,
        )["status"],
    )
    status_b_split = int(
        _classify_topology(
            base=base_b,
            ratio=ratio_b,
            v_piston=case_b["v_piston"],
            energy_split=True,
            root_samples=4000,
        )["status"],
    )
    if status_b_orig != STATUS_THREE_WAVE or status_b_split != STATUS_THREE_WAVE:
        raise RuntimeError("Selected fictitious case is not three-wave in both modes.")

    out_a = assets / "selected_three_wave_case_comparison.png"
    out_b = assets / "fictitious_three_wave_case_comparison.png"

    m_a = _plot_case(
        label=(
            "Selected 3-wave case (map-guided): "
            f"Y0/G={case_a['Y_0']/case_a['G']:.3f}, v_piston={case_a['v_piston']:.0f} cm/s"
        ),
        kwargs=case_a,
        t=2.0e-6,
        out_path=out_a,
    )
    m_b = _plot_case(
        label=(
            "Fictitious 3-wave case (large mode difference): "
            f"Y0/G={case_b['Y_0']/case_b['G']:.3f}, v_piston={case_b['v_piston']:.0f} cm/s"
        ),
        kwargs=case_b,
        t=1.0e-6,
        out_path=out_b,
    )

    summary_path = Path(__file__).resolve().parent / "selected_three_wave_cases.md"
    _build_summary_md(
        summary_path,
        [
            {"label": "Case A", "params": case_a, "metrics": m_a, "image": "../assets/" + out_a.name},
            {"label": "Case B (fictitious)", "params": case_b, "metrics": m_b, "image": "../assets/" + out_b.name},
        ],
    )

    print(f"Saved: {out_a}")
    print(f"Saved: {out_b}")
    print(f"Saved: {summary_path}")


if __name__ == "__main__":
    main()
