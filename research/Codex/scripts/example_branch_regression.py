"""Lightweight regression checks for packaged branch-analysis behavior.

Run:
    python3 research/Codex/scripts/example_branch_regression.py
"""

from __future__ import annotations

import sys
from pathlib import Path

_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from elastoplastic_piston_solver import ElastoplasticPistonSolver
from example_branch_sweep import STATUS_ONE_STRONG_SHOCK, _classify_topology


def _assert_close_baseline_modes() -> None:
    kwargs = dict(
        rho_0=2.79,
        C_0=5.33e5,
        s=1.34,
        Gamma_0=2.0,
        G=2.86e11,
        Y_0=2.6e9,
        e_initial=0.0,
        v_piston=5000.0,
    )
    s_orig = ElastoplasticPistonSolver(**kwargs, energy_split=False)
    s_split = ElastoplasticPistonSolver(**kwargs, energy_split=True)

    rel_us = abs(s_split.U_s - s_orig.U_s) / abs(s_orig.U_s)
    assert rel_us < 1.0e-4, f"Unexpected U_s mismatch in baseline case: {rel_us:.3e}"


def _assert_target_classifies_one_strong_shock() -> None:
    base = dict(
        rho_0=2.79,
        C_0=5.33e5,
        s=1.34,
        Gamma_0=2.0,
        G=2.86e11,
        e_initial=0.0,
    )
    topo_split = _classify_topology(
        base=base,
        ratio=0.15,
        v_piston=80000.0,
        energy_split=True,
        root_samples=12000,
    )
    assert int(topo_split["status"]) == STATUS_ONE_STRONG_SHOCK


def _assert_no_plastic_shock_case() -> None:
    kwargs = dict(
        rho_0=2.79,
        C_0=5.33e5,
        s=1.34,
        Gamma_0=2.0,
        G=2.86e11,
        Y_0=0.15 * 2.86e11,
        e_initial=0.0,
        v_piston=50000.0,  # below v_Y for this setup
    )
    try:
        ElastoplasticPistonSolver(**kwargs, energy_split=False)
    except ValueError as exc:
        msg = str(exc)
        assert "must exceed the elastic-limit particle velocity" in msg
        return
    raise AssertionError("Expected no-plastic-shock error was not raised.")


def main() -> None:
    _assert_close_baseline_modes()
    _assert_target_classifies_one_strong_shock()
    _assert_no_plastic_shock_case()
    print("All branch regression checks passed.")


if __name__ == "__main__":
    main()
