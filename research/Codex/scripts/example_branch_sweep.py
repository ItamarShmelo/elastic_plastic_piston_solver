"""Generate phase maps and report for three-wave vs one-strong-shock regimes.

Run:
    python3 research/Codex/scripts/example_branch_sweep.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
from scipy.optimize import brentq

_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from elastoplastic_piston_solver import _mie_gruneisen_pressure

# Regime labels used in maps/reports.
STATUS_THREE_WAVE = 0
STATUS_ONE_STRONG_SHOCK = 1
STATUS_NO_PLASTIC_SHOCK = 2
STATUS_NO_ROOT = 3

STATUS_LABELS = {
    STATUS_THREE_WAVE: "three_wave",
    STATUS_ONE_STRONG_SHOCK: "one_strong_shock",
    STATUS_NO_PLASTIC_SHOCK: "no_plastic_shock",
    STATUS_NO_ROOT: "no_root",
}

STATUS_COLORS = {
    STATUS_THREE_WAVE: "#2ca02c",
    STATUS_ONE_STRONG_SHOCK: "#d62728",
    STATUS_NO_PLASTIC_SHOCK: "#f2b134",
    STATUS_NO_ROOT: "#7f7f7f",
}

PHASE_GRID_ROOT_SCAN_SAMPLES = 600


def _scan_sign_change_roots(
    func,
    x_min: float,
    x_max: float,
    samples: int,
) -> list[float]:
    xs = np.geomspace(x_min, x_max, samples)
    fs = np.array([func(float(x)) for x in xs], dtype=np.float64)
    roots: list[float] = []
    for i in range(samples - 1):
        f0 = float(fs[i])
        f1 = float(fs[i + 1])
        if not (math.isfinite(f0) and math.isfinite(f1)):
            continue
        if f0 * f1 < 0.0:
            roots.append(float(brentq(func, float(xs[i]), float(xs[i + 1]))))

    roots.sort()
    deduped: list[float] = []
    for root in roots:
        if not deduped or abs(root - deduped[-1]) > 1.0e-6:
            deduped.append(root)
    return deduped


def _classify_topology(
    base: dict[str, float],
    ratio: float,
    v_piston: float,
    energy_split: bool,
    root_samples: int,
) -> dict[str, object]:
    rho_0 = float(base["rho_0"])
    C_0 = float(base["C_0"])
    s = float(base["s"])
    gamma_0 = float(base["Gamma_0"])
    G = float(base["G"])
    e_0 = float(base["e_initial"])
    Y_0 = ratio * G

    rho_Y = rho_0 * math.exp(Y_0 / (2.0 * G))
    two_thirds_Y0 = 2.0 / 3.0 * Y_0
    four_thirds_Y0 = 4.0 / 3.0 * Y_0
    e_el_Y = Y_0**2 / (6.0 * rho_Y * G)

    def eos(e: float, rho: float) -> float:
        return _mie_gruneisen_pressure(e, rho, rho_0, C_0, s, gamma_0)

    e_lo = e_0
    e_hi = e_0 + 10.0 * C_0**2

    try:
        if energy_split:

            def f_yield(e_th: float) -> float:
                P = eos(e_th, rho_Y)
                return (
                    e_th
                    + e_el_Y
                    - e_0
                    - 0.5 / (rho_Y * rho_0) * (P + two_thirds_Y0) * (rho_Y - rho_0)
                )

            e_th_Y = float(brentq(f_yield, e_lo, e_hi))
            P_Y = eos(e_th_Y, rho_Y)
        else:

            def f_yield(e: float) -> float:
                P = eos(e, rho_Y)
                return (
                    e
                    - e_0
                    - 0.5 / (rho_Y * rho_0) * (P + two_thirds_Y0) * (rho_Y - rho_0)
                )

            e_Y = float(brentq(f_yield, e_lo, e_hi))
            P_Y = eos(e_Y, rho_Y)
    except ValueError:
        return {
            "status": STATUS_NO_ROOT,
            "U_se": float("nan"),
            "v_Y": float("nan"),
            "roots_below": [],
            "roots_above": [],
        }

    U_se = math.sqrt(rho_Y * (P_Y + two_thirds_Y0) / (rho_0 * (rho_Y - rho_0)))
    v_Y = (rho_Y - rho_0) / rho_Y * U_se

    if v_piston <= v_Y:
        return {
            "status": STATUS_NO_PLASTIC_SHOCK,
            "U_se": U_se,
            "v_Y": v_Y,
            "roots_below": [],
            "roots_above": [],
        }

    if energy_split:

        def f_shock(U_s: float) -> float:
            P_2 = P_Y + rho_Y * (U_s - v_Y) * (v_piston - v_Y)
            rho_2 = rho_Y * (U_s - v_Y) / (U_s - v_piston)
            e_th_2 = (
                e_th_Y
                + e_el_Y
                + 0.5 / (rho_Y * rho_2) * (P_Y + P_2 + four_thirds_Y0) * (rho_2 - rho_Y)
                - Y_0**2 / (6.0 * rho_2 * G)
            )
            return P_2 - eos(e_th_2, rho_2)
    else:

        def f_shock(U_s: float) -> float:
            P_2 = P_Y + rho_Y * (U_s - v_Y) * (v_piston - v_Y)
            rho_2 = rho_Y * (U_s - v_Y) / (U_s - v_piston)
            e_2 = e_Y + 0.5 / (rho_Y * rho_2) * (P_Y + P_2 + four_thirds_Y0) * (rho_2 - rho_Y)
            return P_2 - eos(e_2, rho_2)

    roots = _scan_sign_change_roots(
        f_shock,
        v_piston * 1.001,
        max(2.0 * U_se, v_piston * 1.2),
        samples=root_samples,
    )
    roots_below = [r for r in roots if r < U_se]
    roots_above = [r for r in roots if r > U_se]

    if roots_above:
        # User-requested interpretation: above-U_se branch => strong-shock regime.
        status = STATUS_ONE_STRONG_SHOCK
    elif roots_below:
        status = STATUS_THREE_WAVE
    else:
        status = STATUS_NO_ROOT

    return {
        "status": status,
        "U_se": U_se,
        "v_Y": v_Y,
        "roots_below": roots_below,
        "roots_above": roots_above,
    }


def _scan_phase(
    base: dict[str, float],
    ratios: np.ndarray,
    v_pistons: np.ndarray,
    energy_split: bool,
) -> np.ndarray:
    grid = np.zeros((len(v_pistons), len(ratios)), dtype=np.int8)
    for i_v, v_piston in enumerate(v_pistons):
        for i_r, ratio in enumerate(ratios):
            topo = _classify_topology(
                base,
                ratio=float(ratio),
                v_piston=float(v_piston),
                energy_split=energy_split,
                root_samples=PHASE_GRID_ROOT_SCAN_SAMPLES,
            )
            grid[i_v, i_r] = int(topo["status"])
    return grid


def _contiguous_intervals(values: np.ndarray, mask: np.ndarray) -> list[tuple[float, float]]:
    intervals: list[tuple[float, float]] = []
    start_idx: int | None = None
    for idx, active in enumerate(mask):
        if active and start_idx is None:
            start_idx = idx
        if start_idx is not None and (not active):
            intervals.append((float(values[start_idx]), float(values[idx - 1])))
            start_idx = None
    if start_idx is not None:
        intervals.append((float(values[start_idx]), float(values[-1])))
    return intervals


def _format_ratio_intervals(intervals: list[tuple[float, float]]) -> str:
    if not intervals:
        return "none"
    return ", ".join(f"[{lo:.4f}, {hi:.4f}]" for lo, hi in intervals)


def _format_speed_intervals(intervals: list[tuple[float, float]]) -> str:
    if not intervals:
        return "none"
    return ", ".join(f"[{lo:.0f}, {hi:.0f}] cm/s" for lo, hi in intervals)


def _status_counts(grid: np.ndarray) -> dict[int, int]:
    return {code: int(np.sum(grid == code)) for code in STATUS_LABELS}


def _plot_phase_maps(
    ratios: np.ndarray,
    v_pistons: np.ndarray,
    grid_orig: np.ndarray,
    grid_split: np.ndarray,
    plot_path: Path,
) -> None:
    cmap = ListedColormap([STATUS_COLORS[i] for i in range(4)])
    norm = BoundaryNorm(np.arange(-0.5, 4.5, 1.0), cmap.N)
    extent = [
        float(ratios[0]),
        float(ratios[-1]),
        float(v_pistons[0] / 1.0e5),
        float(v_pistons[-1] / 1.0e5),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.8), sharey=True)
    for ax, grid, title in (
        (axes[0], grid_orig, "Original Mode"),
        (axes[1], grid_split, "Energy-Split Mode"),
    ):
        ax.imshow(
            grid,
            origin="lower",
            aspect="auto",
            interpolation="nearest",
            extent=extent,
            cmap=cmap,
            norm=norm,
        )
        ax.set_title(title)
        ax.set_xlabel("Y0/G")
        ax.plot(
            0.15,
            0.8,
            marker="*",
            markersize=11,
            markerfacecolor="white",
            markeredgecolor="black",
            linestyle="None",
            label="target case",
        )
        ax.legend(loc="upper right", fontsize=8, frameon=True)

    axes[0].set_ylabel(r"$v_{piston}$ ($10^5$ cm/s)")
    fig.suptitle("Regime Map: Three-Wave vs One-Strong-Shock", y=0.99)

    legend_handles = [
        Patch(color=STATUS_COLORS[STATUS_THREE_WAVE], label="three_wave"),
        Patch(color=STATUS_COLORS[STATUS_ONE_STRONG_SHOCK], label="one_strong_shock"),
        Patch(color=STATUS_COLORS[STATUS_NO_PLASTIC_SHOCK], label="no_plastic_shock"),
        Patch(color=STATUS_COLORS[STATUS_NO_ROOT], label="no_root"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, -0.01),
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.96])
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)


def _build_report(
    ratios: np.ndarray,
    v_pistons: np.ndarray,
    grid_orig: np.ndarray,
    grid_split: np.ndarray,
    target_topology_orig: dict[str, object],
    target_topology_split: dict[str, object],
    report_path: Path,
    plot_path: Path,
) -> None:
    target_ratio = 0.15
    target_v = 8.0e4

    ratio_line = np.linspace(0.01, 1.0, 161)
    ratio_status_split = np.array(
        [
            int(
                _classify_topology(
                    base=dict(rho_0=2.79, C_0=5.33e5, s=1.34, Gamma_0=2.0, G=2.86e11, e_initial=0.0),
                    ratio=float(r),
                    v_piston=target_v,
                    energy_split=True,
                    root_samples=2000,
                )["status"]
            )
            for r in ratio_line
        ],
        dtype=np.int8,
    )

    speed_line = np.linspace(1.0e4, 3.0e5, 121)
    speed_status_split = np.array(
        [
            int(
                _classify_topology(
                    base=dict(rho_0=2.79, C_0=5.33e5, s=1.34, Gamma_0=2.0, G=2.86e11, e_initial=0.0),
                    ratio=target_ratio,
                    v_piston=float(v),
                    energy_split=True,
                    root_samples=2000,
                )["status"]
            )
            for v in speed_line
        ],
        dtype=np.int8,
    )

    ratio_three_wave = _contiguous_intervals(ratio_line, ratio_status_split == STATUS_THREE_WAVE)
    ratio_one_shock = _contiguous_intervals(ratio_line, ratio_status_split == STATUS_ONE_STRONG_SHOCK)
    speed_three_wave = _contiguous_intervals(speed_line, speed_status_split == STATUS_THREE_WAVE)
    speed_one_shock = _contiguous_intervals(speed_line, speed_status_split == STATUS_ONE_STRONG_SHOCK)

    counts_orig = _status_counts(grid_orig)
    counts_split = _status_counts(grid_split)
    total_cells = int(grid_orig.size)

    lines: list[str] = []
    lines.append("# Strong-Shock Regime Investigation")
    lines.append("")
    lines.append("## Setup")
    lines.append("")
    lines.append("- Material family: `rho0=2.79`, `C0=5.33e5`, `s=1.34`, `Gamma0=2.0`, `G=2.86e11` (CGS).")
    lines.append("- Map range: `Y0/G` from `0.01` to `1.0`, `v_piston` from `1e4` to `3e5 cm/s`.")
    lines.append("- Classification rule requested: if any shock root exists above `U_se`, mark as `one_strong_shock`.")
    lines.append("")
    lines.append("## Phase Map")
    lines.append("")
    lines.append(f"![Regime phase map](../assets/{plot_path.name})")
    lines.append("")
    lines.append(
        f"Grid resolution: `{len(ratios)} x {len(v_pistons)}` "
        f"(`Y0/G` x `v_piston`), total `{total_cells}` cells per mode."
    )
    lines.append("")
    lines.append("Status definitions:")
    lines.append("- `three_wave`: roots only below `U_se`.")
    lines.append("- `one_strong_shock`: at least one root above `U_se`.")
    lines.append("- `no_plastic_shock`: `v_piston <= v_Y`.")
    lines.append("- `no_root`: no detected sign-change root in scan range.")
    lines.append("")
    lines.append("Coverage counts:")
    lines.append(
        f"- Original: three_wave={counts_orig[STATUS_THREE_WAVE]}, "
        f"one_strong_shock={counts_orig[STATUS_ONE_STRONG_SHOCK]}, "
        f"no_plastic_shock={counts_orig[STATUS_NO_PLASTIC_SHOCK]}, "
        f"no_root={counts_orig[STATUS_NO_ROOT]}"
    )
    lines.append(
        f"- Split: three_wave={counts_split[STATUS_THREE_WAVE]}, "
        f"one_strong_shock={counts_split[STATUS_ONE_STRONG_SHOCK]}, "
        f"no_plastic_shock={counts_split[STATUS_NO_PLASTIC_SHOCK]}, "
        f"no_root={counts_split[STATUS_NO_ROOT]}"
    )
    lines.append("")
    lines.append("## Split-Mode Windows at Target Slices")
    lines.append("")
    lines.append(
        "At `v_piston = 80000 cm/s`:\n"
        f"- `three_wave` in `Y0/G`: {_format_ratio_intervals(ratio_three_wave)}\n"
        f"- `one_strong_shock` in `Y0/G`: {_format_ratio_intervals(ratio_one_shock)}"
    )
    lines.append("")
    lines.append(
        "At `Y0/G = 0.15`:\n"
        f"- `three_wave` in `v_piston`: {_format_speed_intervals(speed_three_wave)}\n"
        f"- `one_strong_shock` in `v_piston`: {_format_speed_intervals(speed_one_shock)}"
    )
    lines.append("")
    lines.append("## Target Case Root Topology (`Y0/G=0.15`, `v_piston=80000`)")
    lines.append("")
    lines.append(
        f"Original mode: `U_se={float(target_topology_orig['U_se']):.2f}`, "
        f"roots_below={ [round(float(v), 2) for v in target_topology_orig['roots_below']] }, "
        f"roots_above={ [round(float(v), 2) for v in target_topology_orig['roots_above']] }, "
        f"status=`{STATUS_LABELS[int(target_topology_orig['status'])]}`."
    )
    lines.append(
        f"Split mode: `U_se={float(target_topology_split['U_se']):.2f}`, "
        f"roots_below={ [round(float(v), 2) for v in target_topology_split['roots_below']] }, "
        f"roots_above={ [round(float(v), 2) for v in target_topology_split['roots_above']] }, "
        f"status=`{STATUS_LABELS[int(target_topology_split['status'])]}`."
    )
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 research/Codex/scripts/example_branch_sweep.py")
    lines.append("```")
    lines.append("")
    lines.append(f"- Plot: `research/Codex/assets/{plot_path.name}`")
    lines.append(f"- Report: `research/Codex/scripts/{report_path.name}`")
    lines.append("")

    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    assets_dir = Path(__file__).resolve().parent.parent / "assets"
    assets_dir.mkdir(exist_ok=True)

    report_path = Path(__file__).resolve().parent / "branch_transition_sweep.md"
    plot_path = assets_dir / "branch_phase_map.png"

    base = dict(
        rho_0=2.79,
        C_0=5.33e5,
        s=1.34,
        Gamma_0=2.0,
        G=2.86e11,
        e_initial=0.0,
    )

    # User-requested ranges.
    ratios = np.linspace(0.01, 1.0, 61)
    v_pistons = np.linspace(1.0e4, 3.0e5, 61)

    grid_orig = _scan_phase(base, ratios, v_pistons, energy_split=False)
    grid_split = _scan_phase(base, ratios, v_pistons, energy_split=True)

    _plot_phase_maps(ratios, v_pistons, grid_orig, grid_split, plot_path)

    target_topology_orig = _classify_topology(
        base,
        ratio=0.15,
        v_piston=8.0e4,
        energy_split=False,
        root_samples=32000,
    )
    target_topology_split = _classify_topology(
        base,
        ratio=0.15,
        v_piston=8.0e4,
        energy_split=True,
        root_samples=32000,
    )

    _build_report(
        ratios=ratios,
        v_pistons=v_pistons,
        grid_orig=grid_orig,
        grid_split=grid_split,
        target_topology_orig=target_topology_orig,
        target_topology_split=target_topology_split,
        report_path=report_path,
        plot_path=plot_path,
    )

    print(f"Wrote plot: {plot_path}")
    print(f"Wrote report: {report_path}")


if __name__ == "__main__":
    main()
