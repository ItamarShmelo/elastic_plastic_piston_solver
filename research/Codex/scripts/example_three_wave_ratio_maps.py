"""Create zoomed ratio maps in the physical three-wave region.

Ratios are defined as:
    quantity_ratio = quantity_original / quantity_split

Included quantities:
    - Yield region: rho_Y, sigma_Y, v_Y
    - Shocked region: rho_2, sigma_2, v_2
    - Wave speeds: U_se, U_s

Run:
    python3 research/Codex/scripts/example_three_wave_ratio_maps.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from scipy.optimize import brentq

_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from elastoplastic_piston_solver import _mie_gruneisen_pressure
from example_branch_sweep import STATUS_THREE_WAVE, _classify_topology


def _scan_sign_change_roots(func, x_min: float, x_max: float, samples: int) -> list[float]:
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


def _three_wave_state(
    base: dict[str, float],
    ratio: float,
    v_piston: float,
    energy_split: bool,
    root_samples: int,
) -> dict[str, float] | None:
    """Return solved three-wave state, or None if not three-wave."""
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
        return None

    U_se = math.sqrt(rho_Y * (P_Y + two_thirds_Y0) / (rho_0 * (rho_Y - rho_0)))
    v_Y = (rho_Y - rho_0) / rho_Y * U_se

    if v_piston <= v_Y:
        return None

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

    if roots_above or not roots_below:
        return None

    U_s = max(roots_below)
    P_2 = P_Y + rho_Y * (U_s - v_Y) * (v_piston - v_Y)
    rho_2 = rho_Y * (U_s - v_Y) / (U_s - v_piston)

    sigma_Y = -two_thirds_Y0 - P_Y
    sigma_2 = -two_thirds_Y0 - P_2

    return {
        "rho_Y": rho_Y,
        "sigma_Y": sigma_Y,
        "v_Y": v_Y,
        "rho_2": rho_2,
        "sigma_2": sigma_2,
        "v_2": v_piston,
        "U_se": U_se,
        "U_s": U_s,
    }


def main() -> None:
    base = dict(
        rho_0=2.79,
        C_0=5.33e5,
        s=1.34,
        Gamma_0=2.0,
        G=2.86e11,
        e_initial=0.0,
    )

    ratios = np.linspace(0.01, 1.0, 121)
    v_pistons = np.linspace(1.0e4, 3.0e5, 121)

    status_orig = np.zeros((len(v_pistons), len(ratios)), dtype=np.int8)
    status_split = np.zeros((len(v_pistons), len(ratios)), dtype=np.int8)
    for i_v, vp in enumerate(v_pistons):
        for i_r, ratio in enumerate(ratios):
            status_orig[i_v, i_r] = int(
                _classify_topology(
                    base=base,
                    ratio=float(ratio),
                    v_piston=float(vp),
                    energy_split=False,
                    root_samples=500,
                )["status"],
            )
            status_split[i_v, i_r] = int(
                _classify_topology(
                    base=base,
                    ratio=float(ratio),
                    v_piston=float(vp),
                    energy_split=True,
                    root_samples=500,
                )["status"],
            )

    mask_both_three = (status_orig == STATUS_THREE_WAVE) & (status_split == STATUS_THREE_WAVE)
    if not np.any(mask_both_three):
        raise RuntimeError("No overlap where both modes are three-wave.")

    idx_v, idx_r = np.where(mask_both_three)
    pad = 4
    r0 = max(0, int(np.min(idx_r)) - pad)
    r1 = min(len(ratios) - 1, int(np.max(idx_r)) + pad)
    v0 = max(0, int(np.min(idx_v)) - pad)
    v1 = min(len(v_pistons) - 1, int(np.max(idx_v)) + pad)

    ratios_zoom = ratios[r0 : r1 + 1]
    v_zoom = v_pistons[v0 : v1 + 1]

    quantity_keys = [
        "rho_Y",
        "sigma_Y",
        "v_Y",
        "rho_2",
        "sigma_2",
        "v_2",
        "U_se",
        "U_s",
    ]
    quantity_titles = [
        r"$\rho_Y^{orig}/\rho_Y^{split}$",
        r"$\sigma_Y^{orig}/\sigma_Y^{split}$",
        r"$v_Y^{orig}/v_Y^{split}$",
        r"$\rho_2^{orig}/\rho_2^{split}$",
        r"$\sigma_2^{orig}/\sigma_2^{split}$",
        r"$v_2^{orig}/v_2^{split}$",
        r"$U_{se}^{orig}/U_{se}^{split}$",
        r"$U_s^{orig}/U_s^{split}$",
    ]

    ratio_maps = {
        key: np.full((len(v_zoom), len(ratios_zoom)), np.nan, dtype=np.float64)
        for key in quantity_keys
    }

    for i_v, vp in enumerate(v_zoom):
        for i_r, ratio in enumerate(ratios_zoom):
            o = _three_wave_state(base, float(ratio), float(vp), energy_split=False, root_samples=1400)
            s = _three_wave_state(base, float(ratio), float(vp), energy_split=True, root_samples=1400)
            if o is None or s is None:
                continue
            for key in quantity_keys:
                denom = float(s[key])
                if denom == 0.0:
                    continue
                ratio_maps[key][i_v, i_r] = float(o[key]) / denom

    extent = [
        float(ratios_zoom[0]),
        float(ratios_zoom[-1]),
        float(v_zoom[0] / 1.0e5),
        float(v_zoom[-1] / 1.0e5),
    ]

    fig, axes = plt.subplots(2, 4, figsize=(16, 7.6), sharex=True, sharey=True)
    cmap = plt.cm.coolwarm.copy()
    cmap.set_bad(color="#d9d9d9")

    for ax, key, title in zip(axes.ravel(), quantity_keys, quantity_titles, strict=True):
        data = ratio_maps[key]
        finite = data[np.isfinite(data)]
        if finite.size == 0:
            im = ax.imshow(
                np.ma.masked_invalid(data),
                origin="lower",
                aspect="auto",
                interpolation="nearest",
                extent=extent,
                cmap=cmap,
            )
        else:
            dev = float(max(abs(np.min(finite) - 1.0), abs(np.max(finite) - 1.0)))
            if dev < 1.0e-6:
                dev = 1.0e-6
            norm = TwoSlopeNorm(vmin=1.0 - dev, vcenter=1.0, vmax=1.0 + dev)
            im = ax.imshow(
                np.ma.masked_invalid(data),
                origin="lower",
                aspect="auto",
                interpolation="nearest",
                extent=extent,
                cmap=cmap,
                norm=norm,
            )
        ax.set_title(title, fontsize=10)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cbar.ax.tick_params(labelsize=8)

    x_ref = 0.15
    y_ref = 0.8
    if extent[0] <= x_ref <= extent[1] and extent[2] <= y_ref <= extent[3]:
        for ax in axes.ravel():
            ax.plot(
                x_ref,
                y_ref,
                marker="*",
                markersize=9,
                markerfacecolor="white",
                markeredgecolor="black",
                linestyle="None",
            )

    for ax in axes[-1, :]:
        ax.set_xlabel("Y0/G")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$v_{piston}$ ($10^5$ cm/s)")

    fig.suptitle("Three-Wave Region Ratios: Original / Split", fontsize=14, y=0.995)
    fig.tight_layout(rect=[0, 0.02, 1, 0.97])

    assets = Path(__file__).resolve().parent.parent / "assets"
    assets.mkdir(exist_ok=True)
    out = assets / "three_wave_ratio_maps_orig_over_split.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)

    print(f"Saved: {out}")
    print(
        "Zoom bounds: "
        f"Y0/G=[{ratios_zoom[0]:.4f}, {ratios_zoom[-1]:.4f}], "
        f"v_piston=[{v_zoom[0]:.0f}, {v_zoom[-1]:.0f}] cm/s",
    )


if __name__ == "__main__":
    main()
