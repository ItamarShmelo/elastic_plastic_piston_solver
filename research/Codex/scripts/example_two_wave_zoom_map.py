"""Generate zoomed U_se/U_s maps for one-strong-shock and three-wave regions.

Run:
    python3 research/Codex/scripts/example_two_wave_zoom_map.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from example_branch_sweep import (
    STATUS_ONE_STRONG_SHOCK,
    STATUS_THREE_WAVE,
    _classify_topology,
)


def _compute_mode_maps(
    base: dict[str, float],
    ratios: np.ndarray,
    v_pistons: np.ndarray,
    energy_split: bool,
    root_samples: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return status and U_se/U_s maps for one mode."""
    status = np.zeros((len(v_pistons), len(ratios)), dtype=np.int8)
    use_over_us = np.full((len(v_pistons), len(ratios)), np.nan, dtype=np.float64)

    for i_v, vp in enumerate(v_pistons):
        for i_r, ratio in enumerate(ratios):
            topo = _classify_topology(
                base=base,
                ratio=float(ratio),
                v_piston=float(vp),
                energy_split=energy_split,
                root_samples=root_samples,
            )
            s = int(topo["status"])
            status[i_v, i_r] = s

            u_se = float(topo["U_se"])
            roots_below = [float(r) for r in topo["roots_below"]]
            roots_above = [float(r) for r in topo["roots_above"]]

            if s == STATUS_THREE_WAVE and roots_below:
                # Physical 3-wave branch: root nearest U_se from below.
                u_s = max(roots_below)
                use_over_us[i_v, i_r] = u_se / u_s
            elif s == STATUS_ONE_STRONG_SHOCK and roots_above:
                # Strong-shock branch: nearest root above U_se.
                u_s = min(roots_above)
                use_over_us[i_v, i_r] = u_se / u_s

    return status, use_over_us


def _save_zoom_panels(
    *,
    ratios: np.ndarray,
    v_pistons: np.ndarray,
    mask_orig: np.ndarray,
    mask_split: np.ndarray,
    ratio_orig: np.ndarray,
    ratio_split: np.ndarray,
    class_title: str,
    ratio_title: str,
    out_panel: Path,
    out_compare: Path,
) -> None:
    mask = mask_orig | mask_split
    if not np.any(mask):
        raise RuntimeError(f"No cells found for requested zoom: {class_title}")

    idx_v, idx_r = np.where(mask)
    pad = 3
    r0 = max(0, int(np.min(idx_r)) - pad)
    r1 = min(len(ratios) - 1, int(np.max(idx_r)) + pad)
    v0 = max(0, int(np.min(idx_v)) - pad)
    v1 = min(len(v_pistons) - 1, int(np.max(idx_v)) + pad)

    ratios_zoom = ratios[r0 : r1 + 1]
    vp_zoom = v_pistons[v0 : v1 + 1]
    mask_zoom_split = mask_split[v0 : v1 + 1, r0 : r1 + 1]
    ratio_zoom_split = ratio_split[v0 : v1 + 1, r0 : r1 + 1]
    ratio_zoom_orig = ratio_orig[v0 : v1 + 1, r0 : r1 + 1]

    class_data = np.where(mask_zoom_split, 1.0, 0.0)
    class_cmap = ListedColormap(["#d9d9d9", "#d62728"])

    ratio_masked_split = np.ma.masked_invalid(ratio_zoom_split)
    ratio_masked_orig = np.ma.masked_invalid(ratio_zoom_orig)
    ratio_cmap = plt.cm.viridis.copy()
    ratio_cmap.set_bad(color="#d9d9d9")

    extent = [
        float(ratios_zoom[0]),
        float(ratios_zoom[-1]),
        float(vp_zoom[0] / 1.0e5),
        float(vp_zoom[-1] / 1.0e5),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.4), sharey=True)
    axes[0].imshow(
        class_data,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        extent=extent,
        cmap=class_cmap,
        vmin=0.0,
        vmax=1.0,
    )
    axes[0].set_title(class_title)
    axes[0].set_xlabel("Y0/G")
    axes[0].set_ylabel(r"$v_{piston}$ ($10^5$ cm/s)")

    im = axes[1].imshow(
        ratio_masked_split,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        extent=extent,
        cmap=ratio_cmap,
    )
    axes[1].set_title(ratio_title)
    axes[1].set_xlabel("Y0/G")

    x_ref = 0.15
    y_ref = 0.8
    if extent[0] <= x_ref <= extent[1] and extent[2] <= y_ref <= extent[3]:
        for ax in axes:
            ax.plot(
                x_ref,
                y_ref,
                marker="*",
                markersize=11,
                markerfacecolor="white",
                markeredgecolor="black",
                linestyle="None",
            )

    cbar = fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04)
    cbar.set_label(r"$U_{se}/U_s$")
    fig.tight_layout()
    fig.savefig(out_panel, dpi=180)
    plt.close(fig)

    finite_orig = ratio_orig[np.isfinite(ratio_orig)]
    finite_split = ratio_split[np.isfinite(ratio_split)]
    finite_both = np.concatenate([finite_orig, finite_split])
    vmin = float(np.min(finite_both))
    vmax = float(np.max(finite_both))

    fig2, axes2 = plt.subplots(1, 2, figsize=(12.4, 5.4), sharey=True)
    im2 = axes2[0].imshow(
        ratio_masked_orig,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        extent=extent,
        cmap=ratio_cmap,
        vmin=vmin,
        vmax=vmax,
    )
    axes2[0].set_title(r"Original Mode: $U_{se}/U_s$")
    axes2[0].set_xlabel("Y0/G")
    axes2[0].set_ylabel(r"$v_{piston}$ ($10^5$ cm/s)")

    axes2[1].imshow(
        ratio_masked_split,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        extent=extent,
        cmap=ratio_cmap,
        vmin=vmin,
        vmax=vmax,
    )
    axes2[1].set_title(r"Energy-Split Mode: $U_{se}/U_s$")
    axes2[1].set_xlabel("Y0/G")

    if extent[0] <= x_ref <= extent[1] and extent[2] <= y_ref <= extent[3]:
        for ax in axes2:
            ax.plot(
                x_ref,
                y_ref,
                marker="*",
                markersize=11,
                markerfacecolor="white",
                markeredgecolor="black",
                linestyle="None",
            )

    cbar2 = fig2.colorbar(im2, ax=axes2.ravel().tolist(), fraction=0.035, pad=0.03)
    cbar2.set_label(r"$U_{se}/U_s$")
    fig2.tight_layout()
    fig2.savefig(out_compare, dpi=180)
    plt.close(fig2)

    print(f"Saved: {out_panel}")
    print(f"Saved: {out_compare}")
    print(
        "Zoom bounds: "
        f"Y0/G=[{ratios_zoom[0]:.4f}, {ratios_zoom[-1]:.4f}], "
        f"v_piston=[{vp_zoom[0]:.0f}, {vp_zoom[-1]:.0f}] cm/s",
    )


def main() -> None:
    base = dict(
        rho_0=2.79,
        C_0=5.33e5,
        s=1.34,
        Gamma_0=2.0,
        G=2.86e11,
        e_initial=0.0,
    )

    # Requested global ranges
    ratios = np.linspace(0.01, 1.0, 121)
    v_pistons = np.linspace(1.0e4, 3.0e5, 121)

    status_orig, ratio_orig = _compute_mode_maps(
        base=base,
        ratios=ratios,
        v_pistons=v_pistons,
        energy_split=False,
        root_samples=1200,
    )
    status_split, ratio_split = _compute_mode_maps(
        base=base,
        ratios=ratios,
        v_pistons=v_pistons,
        energy_split=True,
        root_samples=1200,
    )

    assets = Path(__file__).resolve().parent.parent / "assets"
    assets.mkdir(exist_ok=True)

    # Earlier "two-wave/strong-shock" zoom.
    _save_zoom_panels(
        ratios=ratios,
        v_pistons=v_pistons,
        mask_orig=status_orig == STATUS_ONE_STRONG_SHOCK,
        mask_split=status_split == STATUS_ONE_STRONG_SHOCK,
        ratio_orig=ratio_orig,
        ratio_split=ratio_split,
        class_title="Zoomed Strong-Shock Area (Split Mode)",
        ratio_title=r"Zoomed $U_{se}/U_s$ in Strong-Shock Area",
        out_panel=assets / "two_wave_zoom_use_over_us.png",
        out_compare=assets / "two_wave_zoom_use_over_us_compare_modes.png",
    )

    # Corrected "three-wave" zoom.
    _save_zoom_panels(
        ratios=ratios,
        v_pistons=v_pistons,
        mask_orig=status_orig == STATUS_THREE_WAVE,
        mask_split=status_split == STATUS_THREE_WAVE,
        ratio_orig=ratio_orig,
        ratio_split=ratio_split,
        class_title="Zoomed Physical 3-Wave Area (Split Mode)",
        ratio_title=r"Zoomed $U_{se}/U_s$ in 3-Wave Area",
        out_panel=assets / "three_wave_zoom_use_over_us.png",
        out_compare=assets / "three_wave_zoom_use_over_us_compare_modes.png",
    )


if __name__ == "__main__":
    main()
