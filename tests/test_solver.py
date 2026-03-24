"""Tests for the elastoplastic piston solver.

Covers baseline regression, regime classification, EOS contract validation,
edge cases, and output shape/key checks.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from elastoplastic_piston_solver import (
    ElastoplasticPistonSolver,
    WaveStructure,
)

# ---------------------------------------------------------------------------
# Material parameter sets
# ---------------------------------------------------------------------------

ALUMINUM = dict(
    rho_0=2.79, C_0=5.33e5, s=1.34, Gamma_0=2.0,
    G=2.86e11, Y_0=2.6e9, e_initial=0.0,
)

COPPER = dict(
    rho_0=8.93, C_0=3.94e5, s=1.489, Gamma_0=2.02,
    G=4.77e11, Y_0=1.2e9, e_initial=0.0,
)


# ---------------------------------------------------------------------------
# Baseline regression tests
# ---------------------------------------------------------------------------

class TestBaselineRegression:
    """Current aluminum and copper two-wave states must match numerically."""

    def test_aluminum_no_split(self):
        s = ElastoplasticPistonSolver(**ALUMINUM, v_piston=5000.0)
        assert s.wave_structure == WaveStructure.PLASTIC_SHOCK_AND_ELASTIC_WAVE
        assert s.rho_Y == pytest.approx(2.8027106842157106, rel=1e-10)
        assert s.e_Y == pytest.approx(4372546.152088927, rel=1e-10)
        assert s.P_Y == pytest.approx(3646612968.7048893, rel=1e-10)
        assert s.U_se == pytest.approx(652065.9743698481, rel=1e-10)
        assert s.v_Y == pytest.approx(2957.210223196493, rel=1e-10)
        assert s.U_s == pytest.approx(543899.5692549741, rel=1e-6)
        assert s.rho_2 == pytest.approx(2.813334832126923, rel=1e-6)
        assert s.P_2 == pytest.approx(6743696618.642619, rel=1e-6)
        assert s.e_2 == pytest.approx(13707945.33319463, rel=1e-6)
        assert s.equivalent_plastic_strain == pytest.approx(
            0.003783502307422039, rel=1e-6)

    def test_aluminum_split(self):
        s = ElastoplasticPistonSolver(**ALUMINUM, v_piston=5000.0,
                                     energy_split=True)
        assert s.wave_structure == WaveStructure.PLASTIC_SHOCK_AND_ELASTIC_WAVE
        assert s.rho_Y == pytest.approx(2.8027106842157106, rel=1e-10)
        assert s.e_Y == pytest.approx(4366103.610653991, rel=1e-10)
        assert s.P_Y == pytest.approx(3638686116.531571, rel=1e-10)
        assert s.U_se == pytest.approx(651585.4179062786, rel=1e-10)
        assert s.v_Y == pytest.approx(2955.030832547772, rel=1e-10)
        assert s.U_s == pytest.approx(543892.4270612809, rel=1e-6)
        assert s.rho_2 == pytest.approx(2.8133463076661367, rel=1e-6)
        assert s.P_2 == pytest.approx(6739045507.282562, rel=1e-6)
        assert s.e_2 == pytest.approx(13703065.334690029, rel=1e-6)
        assert s.e_th_Y == pytest.approx(2958405.8667642735, rel=1e-10)
        assert s.e_el_Y == pytest.approx(1407697.7438897172, rel=1e-10)
        assert s.e_th_2 == pytest.approx(12295367.590800311, rel=1e-6)
        assert s.e_el_2 == pytest.approx(1407697.7438897172, rel=1e-6)
        assert s.equivalent_plastic_strain == pytest.approx(
            0.003787581280062725, rel=1e-6)

    def test_copper_no_split(self):
        s = ElastoplasticPistonSolver(**COPPER, v_piston=3000.0)
        assert s.wave_structure == WaveStructure.PLASTIC_SHOCK_AND_ELASTIC_WAVE
        assert s.rho_Y == pytest.approx(8.941239771959504, rel=1e-10)
        assert s.e_Y == pytest.approx(179493.9623954217, rel=1e-10)
        assert s.P_Y == pytest.approx(1750184140.9820726, rel=1e-10)
        assert s.U_se == pytest.approx(476628.6690728114, rel=1e-10)
        assert s.v_Y == pytest.approx(599.1560103936564, rel=1e-10)
        assert s.U_s == pytest.approx(399376.5419241198, rel=1e-6)
        assert s.rho_2 == pytest.approx(8.995396664448483, rel=1e-6)
        assert s.P_2 == pytest.approx(10310547575.542755, rel=1e-6)
        assert s.e_2 == pytest.approx(4778665.872251198, rel=1e-6)
        assert s.equivalent_plastic_strain == pytest.approx(
            0.006038708233245766, rel=1e-6)


# ---------------------------------------------------------------------------
# Elastic-only regime tests
# ---------------------------------------------------------------------------

class TestElasticOnlyRegime:

    @pytest.mark.parametrize("v_piston", [100.0, 500.0, 1000.0])
    def test_sub_yield_returns_elastic_wave(self, v_piston):
        s = ElastoplasticPistonSolver(**ALUMINUM, v_piston=v_piston)
        assert s.wave_structure == WaveStructure.ELASTIC_WAVE
        assert s.v_Y == v_piston

    @pytest.mark.parametrize("v_piston", [100.0, 500.0, 1000.0])
    def test_sub_yield_state(self, v_piston):
        """Elastic solve produces sub-yield compression and stress."""
        s = ElastoplasticPistonSolver(**ALUMINUM, v_piston=v_piston)
        s_ref = ElastoplasticPistonSolver(**ALUMINUM, v_piston=5000.0)
        rho_Y_yield = s_ref.rho_Y
        assert s.rho_Y < rho_Y_yield
        assert abs(s.Sx_precursor) < 2.0 / 3.0 * ALUMINUM["Y_0"]

    def test_near_threshold(self):
        s = ElastoplasticPistonSolver(**ALUMINUM, v_piston=2950.0)
        assert s.wave_structure == WaveStructure.ELASTIC_WAVE
        assert s.v_Y == s.v_piston

    def test_continuity_at_yield_boundary(self):
        """At v_piston just below v_Y, elastic state approximates yield."""
        s_ref = ElastoplasticPistonSolver(**ALUMINUM, v_piston=5000.0)
        v_Y_ref = s_ref.v_Y
        rho_Y_ref = s_ref.rho_Y
        P_Y_ref = s_ref.P_Y
        U_se_ref = s_ref.U_se

        s_near = ElastoplasticPistonSolver(
            **ALUMINUM, v_piston=v_Y_ref * (1.0 - 1e-3),
        )
        assert s_near.wave_structure == WaveStructure.ELASTIC_WAVE
        assert s_near.rho_Y == pytest.approx(rho_Y_ref, rel=1e-2)
        assert s_near.P_Y == pytest.approx(P_Y_ref, rel=2e-3)
        assert s_near.U_se == pytest.approx(U_se_ref, rel=1e-2)
        assert s_near.Sx_precursor == pytest.approx(
            -2.0 / 3.0 * ALUMINUM["Y_0"], rel=2e-3,
        )

    def test_solve_output_elastic(self):
        s = ElastoplasticPistonSolver(**ALUMINUM, v_piston=100.0)
        r = s.solve(1e-6, np.linspace(0, 1, 50))
        assert math.isnan(r["shock_location"])
        assert math.isfinite(r["elastic_precursor_location"])
        assert r["wave_structure"] == WaveStructure.ELASTIC_WAVE

    def test_state_2_is_none(self):
        s = ElastoplasticPistonSolver(**ALUMINUM, v_piston=100.0)
        assert s.U_s is None
        assert s.rho_2 is None
        assert s.P_2 is None
        assert s.e_2 is None
        assert s.equivalent_plastic_strain is None

    def test_v_piston_zero_rejected(self):
        with pytest.raises(ValueError, match="v_piston must be positive"):
            ElastoplasticPistonSolver(**ALUMINUM, v_piston=0.0)

    def test_v_piston_negative_rejected(self):
        with pytest.raises(ValueError, match="v_piston must be positive"):
            ElastoplasticPistonSolver(**ALUMINUM, v_piston=-100.0)

    @pytest.mark.parametrize("split", [False, True])
    def test_energy_split_elastic(self, split):
        s = ElastoplasticPistonSolver(**ALUMINUM, v_piston=100.0,
                                     energy_split=split)
        assert s.wave_structure == WaveStructure.ELASTIC_WAVE
        if split:
            a = math.log(s.rho_Y / s.rho_0)
            e_el_expected = (4.0 * s.G) / (3.0 * s.rho_0) * (
                1.0 - (1.0 + a) * math.exp(-a))
            assert s.e_el_Y == pytest.approx(e_el_expected, rel=1e-10)
            assert s.e_th_Y + s.e_el_Y == pytest.approx(s.e_Y, rel=1e-10)

    def test_ideal_gas_elastic(self):
        s = ElastoplasticPistonSolver(
            rho_0=2.79, gamma_ideal_gas=3.0,
            G=2.86e11, Y_0=2.6e9, e_initial=0.0, v_piston=100.0,
        )
        assert s.wave_structure == WaveStructure.ELASTIC_WAVE
        assert s.v_Y == 100.0
        assert math.isfinite(s.rho_Y)
        assert math.isfinite(s.P_Y)
        assert math.isfinite(s.U_se)
        assert abs(s.Sx_precursor) < 2.0 / 3.0 * 2.6e9


# ---------------------------------------------------------------------------
# Two-wave regime tests
# ---------------------------------------------------------------------------

class TestTwoWaveRegime:

    def test_standard_aluminum(self):
        s = ElastoplasticPistonSolver(**ALUMINUM, v_piston=5000.0)
        assert s.wave_structure == WaveStructure.PLASTIC_SHOCK_AND_ELASTIC_WAVE
        assert s.U_s < s.U_se

    def test_standard_copper(self):
        s = ElastoplasticPistonSolver(**COPPER, v_piston=3000.0)
        assert s.wave_structure == WaveStructure.PLASTIC_SHOCK_AND_ELASTIC_WAVE
        assert s.U_s < s.U_se

    @pytest.mark.parametrize("split", [False, True])
    def test_monotonic_compression(self, split):
        s = ElastoplasticPistonSolver(**ALUMINUM, v_piston=5000.0,
                                     energy_split=split)
        assert s.rho_2 > s.rho_Y > s.rho_0

    @pytest.mark.parametrize("split", [False, True])
    def test_equivalent_plastic_strain(self, split):
        s = ElastoplasticPistonSolver(**ALUMINUM, v_piston=5000.0,
                                     energy_split=split)
        assert s.equivalent_plastic_strain > 0
        assert s.equivalent_plastic_strain == pytest.approx(
            math.log(s.rho_2 / s.rho_Y), rel=1e-12)

    def test_solve_output_two_wave(self):
        s = ElastoplasticPistonSolver(**ALUMINUM, v_piston=5000.0)
        r = s.solve(1e-6, np.linspace(0, 1, 100))
        assert math.isfinite(r["shock_location"])
        assert math.isfinite(r["elastic_precursor_location"])
        assert r["shock_location"] < r["elastic_precursor_location"]
        assert r["wave_structure"] == WaveStructure.PLASTIC_SHOCK_AND_ELASTIC_WAVE


# ---------------------------------------------------------------------------
# One plastic wave (plastic-only) regime tests
# ---------------------------------------------------------------------------

class TestPlasticWaveRegime:

    def test_high_velocity_high_strength_split(self):
        """The known problematic case that triggered this refactor."""
        G = 2.86e11
        Y_0 = 0.15 * G
        s = ElastoplasticPistonSolver(
            rho_0=2.79, C_0=5.33e5, s=1.34, Gamma_0=2.0,
            G=G, Y_0=Y_0, e_initial=0.0, v_piston=80000.0,
            energy_split=True,
        )
        assert s.wave_structure == WaveStructure.PLASTIC_WAVE
        assert math.isfinite(s.U_s)
        assert s.rho_2 > s.rho_0

    def test_high_velocity_high_strength_no_split(self):
        """In no-split mode this point still has 2 roots below U_se
        (the energy-split boundary shift only affects split mode)."""
        G = 2.86e11
        Y_0 = 0.15 * G
        s = ElastoplasticPistonSolver(
            rho_0=2.79, C_0=5.33e5, s=1.34, Gamma_0=2.0,
            G=G, Y_0=Y_0, e_initial=0.0, v_piston=80000.0,
            energy_split=False,
        )
        assert s.wave_structure == WaveStructure.PLASTIC_SHOCK_AND_ELASTIC_WAVE
        assert math.isfinite(s.U_s)
        assert s.rho_2 > s.rho_0

    def test_solve_output_plastic_wave(self):
        G = 2.86e11
        Y_0 = 0.15 * G
        s = ElastoplasticPistonSolver(
            rho_0=2.79, C_0=5.33e5, s=1.34, Gamma_0=2.0,
            G=G, Y_0=Y_0, e_initial=0.0, v_piston=80000.0,
            energy_split=True,
        )
        r = s.solve(1e-6, np.linspace(0, 1, 100))
        assert math.isfinite(r["shock_location"])
        assert math.isnan(r["elastic_precursor_location"])
        assert r["wave_structure"] == WaveStructure.PLASTIC_WAVE
        assert "e_thermal" in r
        assert "e_elastic" in r

    @pytest.mark.parametrize("v_piston", [90000.0, 120000.0, 200000.0])
    def test_various_high_velocities(self, v_piston):
        G = 2.86e11
        Y_0 = 0.15 * G
        s = ElastoplasticPistonSolver(
            rho_0=2.79, C_0=5.33e5, s=1.34, Gamma_0=2.0,
            G=G, Y_0=Y_0, e_initial=0.0, v_piston=v_piston,
            energy_split=True,
        )
        assert s.wave_structure == WaveStructure.PLASTIC_WAVE
        assert math.isfinite(s.U_s)
        assert math.isfinite(s.P_2)
        assert s.rho_2 > s.rho_0


# ---------------------------------------------------------------------------
# EOS contract tests
# ---------------------------------------------------------------------------

class TestEOSContract:

    def test_accepts_complete_mg(self):
        s = ElastoplasticPistonSolver(
            rho_0=2.79, C_0=5.33e5, s=1.34, Gamma_0=2.0,
            G=2.86e11, Y_0=2.6e9, e_initial=0.0, v_piston=5000.0,
        )
        assert s.C_0 == 5.33e5

    def test_accepts_gamma_ideal_gas(self):
        s = ElastoplasticPistonSolver(
            rho_0=2.79, gamma_ideal_gas=3.0,
            G=2.86e11, Y_0=2.6e9, e_initial=0.0, v_piston=5000.0,
        )
        assert s.C_0 == 0.0
        assert s.s == 0.0
        assert s.Gamma_0 == 2.0

    def test_rejects_mixed(self):
        with pytest.raises(ValueError, match="mutually exclusive"):
            ElastoplasticPistonSolver(
                rho_0=2.79, C_0=5.33e5, gamma_ideal_gas=3.0,
                G=2.86e11, Y_0=2.6e9, e_initial=0.0, v_piston=5000.0,
            )

    def test_rejects_partial_mg(self):
        with pytest.raises(ValueError, match="Partial"):
            ElastoplasticPistonSolver(
                rho_0=2.79, C_0=5.33e5,
                G=2.86e11, Y_0=2.6e9, e_initial=0.0, v_piston=5000.0,
            )

    def test_rejects_no_eos(self):
        with pytest.raises(ValueError, match="No EOS"):
            ElastoplasticPistonSolver(
                rho_0=2.79,
                G=2.86e11, Y_0=2.6e9, e_initial=0.0, v_piston=5000.0,
            )


# ---------------------------------------------------------------------------
# Ideal-gas behaviour tests
# ---------------------------------------------------------------------------

class TestIdealGas:

    @pytest.mark.parametrize("v_piston", [5000.0, 20000.0])
    def test_ideal_gas_finite_state(self, v_piston):
        s = ElastoplasticPistonSolver(
            rho_0=2.79, gamma_ideal_gas=3.0,
            G=2.86e11, Y_0=2.6e9, e_initial=0.0, v_piston=v_piston,
        )
        assert s.wave_structure in (
            WaveStructure.PLASTIC_SHOCK_AND_ELASTIC_WAVE,
            WaveStructure.PLASTIC_WAVE,
        )
        assert math.isfinite(s.U_s)
        assert math.isfinite(s.P_2)
        assert s.rho_2 > s.rho_0

    def test_ideal_gas_split(self):
        s = ElastoplasticPistonSolver(
            rho_0=2.79, gamma_ideal_gas=3.0,
            G=2.86e11, Y_0=2.6e9, e_initial=0.0, v_piston=5000.0,
            energy_split=True,
        )
        assert math.isfinite(s.U_s)
        assert math.isfinite(s.e_2)

    def test_ideal_gas_two_wave(self):
        """Ideal gas with moderate piston velocity gives two-wave regime."""
        s = ElastoplasticPistonSolver(
            rho_0=2.0, gamma_ideal_gas=3.0,
            G=1e8, Y_0=1e7, e_initial=1e6, v_piston=1000.0,
        )
        assert s.wave_structure == WaveStructure.PLASTIC_SHOCK_AND_ELASTIC_WAVE
        assert s.U_s < s.U_se
        assert math.isfinite(s.U_s)
        assert math.isfinite(s.P_2)
        assert s.rho_2 > s.rho_0

    def test_ideal_gas_plastic_wave(self):
        """Ideal gas with high piston velocity gives plastic wave regime."""
        s = ElastoplasticPistonSolver(
            rho_0=2.0, gamma_ideal_gas=3.0,
            G=1e8, Y_0=1e7, e_initial=1e6, v_piston=10000.0,
        )
        assert s.wave_structure == WaveStructure.PLASTIC_WAVE
        assert math.isfinite(s.U_s)
        assert math.isfinite(s.P_2)
        assert s.rho_2 > s.rho_0


# ---------------------------------------------------------------------------
# Mode coverage (energy_split cross-product)
# ---------------------------------------------------------------------------

class TestModeCoverage:

    @pytest.mark.parametrize("split", [False, True])
    def test_two_wave_modes(self, split):
        s = ElastoplasticPistonSolver(**ALUMINUM, v_piston=5000.0,
                                     energy_split=split)
        assert s.wave_structure == WaveStructure.PLASTIC_SHOCK_AND_ELASTIC_WAVE

    @pytest.mark.parametrize("split", [False, True])
    def test_elastic_modes(self, split):
        s = ElastoplasticPistonSolver(**ALUMINUM, v_piston=100.0,
                                     energy_split=split)
        assert s.wave_structure == WaveStructure.ELASTIC_WAVE

    def test_plastic_wave_mode_split(self):
        G = 2.86e11
        Y_0 = 0.15 * G
        s = ElastoplasticPistonSolver(
            rho_0=2.79, C_0=5.33e5, s=1.34, Gamma_0=2.0,
            G=G, Y_0=Y_0, e_initial=0.0, v_piston=80000.0,
            energy_split=True,
        )
        assert s.wave_structure == WaveStructure.PLASTIC_WAVE

    def test_plastic_wave_mode_no_split_higher_vp(self):
        """Need higher v_piston for no-split to enter one-wave regime."""
        G = 2.86e11
        Y_0 = 0.15 * G
        s = ElastoplasticPistonSolver(
            rho_0=2.79, C_0=5.33e5, s=1.34, Gamma_0=2.0,
            G=G, Y_0=Y_0, e_initial=0.0, v_piston=120000.0,
            energy_split=False,
        )
        assert s.wave_structure == WaveStructure.PLASTIC_WAVE


# ---------------------------------------------------------------------------
# Output shape and key checks
# ---------------------------------------------------------------------------

class TestOutputContract:

    @pytest.fixture(params=[
        ("two_wave", dict(**ALUMINUM, v_piston=5000.0)),
        ("elastic", dict(**ALUMINUM, v_piston=100.0)),
        ("plastic", dict(rho_0=2.79, C_0=5.33e5, s=1.34, Gamma_0=2.0,
                         G=2.86e11, Y_0=0.15*2.86e11, e_initial=0.0,
                         v_piston=80000.0, energy_split=True)),
    ])
    def solver_and_result(self, request):
        label, kwargs = request.param
        s = ElastoplasticPistonSolver(**kwargs)
        r = s.solve(1e-6, np.linspace(0, 1, 50))
        return label, s, r

    def test_required_keys(self, solver_and_result):
        _, _, r = solver_and_result
        for key in ("density", "pressure", "velocity", "energy", "Sx",
                     "stress", "equivalent_plastic_strain",
                     "shock_location", "elastic_precursor_location",
                     "piston_location", "wave_structure"):
            assert key in r

    def test_array_shapes(self, solver_and_result):
        _, _, r = solver_and_result
        for key in ("density", "pressure", "velocity", "energy", "Sx",
                     "stress", "equivalent_plastic_strain"):
            assert r[key].shape == (50,)

    def test_wave_structure_is_enum(self, solver_and_result):
        _, _, r = solver_and_result
        assert isinstance(r["wave_structure"], WaveStructure)

    def test_finite_field_values(self, solver_and_result):
        """Values at and ahead of the piston must be finite; behind is NaN."""
        _, _, r = solver_and_result
        x = np.linspace(0, 1, 50)
        ahead = x >= r["piston_location"]
        for key in ("density", "pressure", "velocity", "energy",
                     "equivalent_plastic_strain"):
            assert np.all(np.isfinite(r[key][ahead]))
        behind = x < r["piston_location"]
        if np.any(behind):
            for key in ("density", "pressure", "velocity", "energy",
                         "equivalent_plastic_strain"):
                assert np.all(np.isnan(r[key][behind]))

    def test_elastic_regime_zero_plastic_strain(self, solver_and_result):
        label, _, r = solver_and_result
        if label != "elastic":
            pytest.skip("only applies to elastic regime")
        x = np.linspace(0, 1, 50)
        ahead = x >= r["piston_location"]
        assert np.all(r["equivalent_plastic_strain"][ahead] == 0.0)


# ---------------------------------------------------------------------------
# Numerical sanity
# ---------------------------------------------------------------------------

class TestNumericalSanity:

    @pytest.mark.parametrize("kwargs", [
        dict(**ALUMINUM, v_piston=5000.0),
        dict(**ALUMINUM, v_piston=5000.0, energy_split=True),
        dict(**COPPER, v_piston=3000.0),
        dict(rho_0=2.79, C_0=5.33e5, s=1.34, Gamma_0=2.0,
             G=2.86e11, Y_0=0.15*2.86e11, e_initial=0.0, v_piston=80000.0),
        dict(rho_0=2.79, gamma_ideal_gas=3.0,
             G=2.86e11, Y_0=2.6e9, e_initial=0.0, v_piston=5000.0),
    ])
    def test_all_solved_quantities_finite(self, kwargs):
        s = ElastoplasticPistonSolver(**kwargs)
        for name in ("rho_Y", "e_Y", "P_Y", "U_se", "v_Y"):
            assert math.isfinite(getattr(s, name))
        if s.wave_structure != WaveStructure.ELASTIC_WAVE:
            for name in ("U_s", "rho_2", "P_2", "e_2",
                         "equivalent_plastic_strain"):
                assert math.isfinite(getattr(s, name))
