"""Analytic solver for the 1D elastoplastic piston test problem.

Solves the problem of a piston driving into an elastoplastic medium with a
Mie-Gruneisen equation of state.  Depending on the input parameters, the
solution consists of:

- An elastic wave only (piston velocity at or below yield threshold),
- A plastic shock preceded by an elastic precursor (two-wave regime), or
- A single plastic shock with no separate elastic precursor (one-wave regime).

Reference
---------
See the derivation in ``README.md`` for the full notation and equations.
"""

from __future__ import annotations

import enum
import math

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import brentq


# ---------------------------------------------------------------------------
# Wave structure enum
# ---------------------------------------------------------------------------

class WaveStructure(enum.Enum):
    """Classification of the wave structure produced by the piston."""

    ELASTIC_WAVE = "elastic wave"
    PLASTIC_SHOCK_AND_ELASTIC_WAVE = "plastic shock and elastic wave"
    PLASTIC_WAVE = "plastic wave"


# ---------------------------------------------------------------------------
# Mie-Gruneisen EOS helper functions
# ---------------------------------------------------------------------------

def _hugoniot_pressure(
    rho: float,
    rho_0: float,
    C_0: float,
    s: float,
) -> float:
    """Hugoniot reference pressure P_H(rho).

    Parameters
    ----------
    rho : float
        Current density.
    rho_0 : float
        Reference density.
    C_0 : float
        Bulk sound speed at reference state.
    s : float
        Hugoniot slope coefficient.

    Returns
    -------
    float
        Hugoniot pressure P_H.
    """
    mu: float = 1.0 - rho_0 / rho
    if rho >= rho_0:
        return rho_0 * C_0**2 * mu / (1.0 - s * mu) ** 2
    else:
        return rho_0 * C_0**2 * mu


def _hugoniot_energy(
    rho: float,
    rho_0: float,
    C_0: float,
    s: float,
) -> float:
    """Hugoniot reference energy e_H(rho).

    e_H = P_H * mu / (2 * rho_0)

    Parameters
    ----------
    rho, rho_0, C_0, s :
        Same as :func:`_hugoniot_pressure`.

    Returns
    -------
    float
        Hugoniot energy e_H.
    """
    mu: float = 1.0 - rho_0 / rho
    P_H: float = _hugoniot_pressure(rho, rho_0, C_0, s)
    return 0.5 * P_H * mu / rho_0


def _mie_gruneisen_pressure(
    e: float,
    rho: float,
    rho_0: float,
    C_0: float,
    s: float,
    Gamma_0: float,
) -> float:
    """Mie-Gruneisen equation of state pressure.

    P(e, rho) = P_H(rho) + Gamma_0 * rho * (e - e_H(rho))

    Parameters
    ----------
    e : float
        Specific internal energy.
    rho : float
        Current density.
    rho_0 : float
        Reference density.
    C_0 : float
        Bulk sound speed at reference state.
    s : float
        Hugoniot slope coefficient.
    Gamma_0 : float
        Gruneisen parameter.

    Returns
    -------
    float
        Pressure.
    """
    P_H: float = _hugoniot_pressure(rho, rho_0, C_0, s)
    e_H: float = _hugoniot_energy(rho, rho_0, C_0, s)
    return P_H + Gamma_0 * rho * (e - e_H)


# ---------------------------------------------------------------------------
# Root-scanning utility
# ---------------------------------------------------------------------------

def _scan_sign_change_roots(
    func,
    x_min: float,
    x_max: float,
    samples: int = 600,
) -> list[float]:
    """Find all roots of *func* in ``[x_min, x_max]`` via sign-change scan.

    Uses geometrically spaced probes, refines each sign-change bracket with
    Brent's method, and deduplicates roots closer than ``1e-6``.
    """
    xs = np.geomspace(x_min, x_max, samples)
    fs = np.array([func(float(x)) for x in xs], dtype=np.float64)
    roots: list[float] = []
    for i in range(samples - 1):
        f0, f1 = float(fs[i]), float(fs[i + 1])
        if math.isfinite(f0) and math.isfinite(f1) and f0 * f1 < 0.0:
            try:
                roots.append(float(brentq(func, float(xs[i]), float(xs[i + 1]))))
            except ValueError:
                pass
    roots.sort()
    deduped: list[float] = []
    for root in roots:
        if not deduped or abs(root - deduped[-1]) > 1.0e-6:
            deduped.append(root)
    return deduped


# ---------------------------------------------------------------------------
# Solver class
# ---------------------------------------------------------------------------

class ElastoplasticPistonSolver:
    """Analytic solver for the 1D elastoplastic piston problem.

    The solver computes the piecewise-constant solution produced when a piston
    drives into an elastoplastic medium at constant velocity.  The resulting
    wave structure is classified as one of:

    * ``WaveStructure.ELASTIC_WAVE`` -- piston velocity at or below the
      elastic-limit particle velocity; only an elastic precursor propagates.
    * ``WaveStructure.PLASTIC_SHOCK_AND_ELASTIC_WAVE`` -- a plastic shock
      trailed by an elastic precursor (the standard two-wave regime).
    * ``WaveStructure.PLASTIC_WAVE`` -- a single plastic shock with no
      separate elastic precursor (strong-shock / one-wave regime).  In this
      regime the Rankine-Hugoniot jump is taken directly from the initial
      state to the shocked state.

    Parameters
    ----------
    rho_0 : float
        Reference density (always required).
    C_0 : float, optional
        Bulk sound speed at reference state (Mie-Gruneisen EOS).
    s : float, optional
        Hugoniot slope coefficient (Mie-Gruneisen EOS).
    Gamma_0 : float, optional
        Gruneisen parameter (Mie-Gruneisen EOS).
    G : float
        Shear modulus.
    Y_0 : float
        Yield strength.
    e_initial : float
        Initial specific internal energy.
    v_piston : float
        Piston velocity (positive, in the +x direction).
    energy_split : bool
        If ``True``, split the internal energy into thermal and elastic
        components and feed only the thermal part to the EOS.  The elastic
        energy is computed from the work integral of deviatoric power.
        Default is ``False``.
    gamma_ideal_gas : float, optional
        Ideal-gas adiabatic index.  Mutually exclusive with ``C_0``, ``s``,
        ``Gamma_0``.  Internally mapped to Mie-Gruneisen constants
        ``C_0 = 0``, ``s = 0``, ``Gamma_0 = gamma_ideal_gas - 1``.

    Raises
    ------
    ValueError
        If the EOS specification is invalid (partial MG tuple, mixed
        ideal-gas + MG, or nothing provided).
    """

    def __init__(
        self,
        rho_0: float,
        C_0: float | None = None,
        s: float | None = None,
        Gamma_0: float | None = None,
        *,
        G: float,
        Y_0: float,
        e_initial: float,
        v_piston: float,
        energy_split: bool = False,
        gamma_ideal_gas: float | None = None,
    ) -> None:
        # -- EOS validation and configuration --------------------------------
        mg_provided = (C_0 is not None, s is not None, Gamma_0 is not None)
        if gamma_ideal_gas is not None:
            if any(mg_provided):
                raise ValueError(
                    "gamma_ideal_gas is mutually exclusive with C_0, s, "
                    "Gamma_0.  Provide either gamma_ideal_gas alone or all "
                    "three Mie-Gruneisen parameters."
                )
            C_0 = 0.0
            s = 0.0
            Gamma_0 = gamma_ideal_gas - 1.0
        elif any(mg_provided):
            if not all(mg_provided):
                raise ValueError(
                    "Partial Mie-Gruneisen EOS specification: C_0, s, and "
                    "Gamma_0 must all be provided together."
                )
        else:
            raise ValueError(
                "No EOS specified.  Provide either gamma_ideal_gas or all "
                "of (C_0, s, Gamma_0)."
            )

        if v_piston <= 0:
            raise ValueError(
                f"v_piston must be positive, got {v_piston}."
            )

        self.rho_0: float = rho_0
        self.C_0: float = C_0
        self.s: float = s
        self.Gamma_0: float = Gamma_0

        self.G: float = G
        self.Y_0: float = Y_0

        self.e_initial: float = e_initial
        self.v_piston: float = v_piston

        self.energy_split: bool = energy_split

        # -- Solved state variables (populated by _solve_wave_structure) -----
        self.wave_structure: WaveStructure | None = None

        self.rho_Y: float | None = None
        self.e_Y: float | None = None
        self.P_Y: float | None = None
        self.U_se: float | None = None
        self.v_Y: float | None = None
        self.U_s: float | None = None
        self.rho_2: float | None = None
        self.P_2: float | None = None
        self.e_2: float | None = None

        self.e_th_Y: float | None = None
        self.e_el_Y: float | None = None
        self.e_th_2: float | None = None
        self.e_el_2: float | None = None

        self.Sx_precursor: float | None = None

        self.P_initial: float = self._eos(e_initial, rho_0)

        self._solve_wave_structure()

    # -- convenience wrappers ----------------------------------------------

    def _eos(self, e: float, rho: float) -> float:
        """Evaluate the Mie-Gruneisen EOS with stored parameters."""
        return _mie_gruneisen_pressure(
            e, rho, self.rho_0, self.C_0, self.s, self.Gamma_0,
        )

    def _elastic_energy_work_integral(self, rho: float) -> float:
        """Elastic energy from the integral of deviatoric power S:D^el.

        Computes  e_el(rho) = (4G)/(3 rho_0) * [1 - (1 + a) exp(-a)]

        where a = ln(rho / rho_0).  This is the integral of
        (-S_x) d(rho) / rho^2 from rho_0 to rho, with the hypoelastic
        constitutive law S_x = (4/3) G ln(rho_0 / rho).

        For small a (typical metals, Y0/(2G) ~ 1e-3), this reduces to
        the state function S:S/(4 rho G) to leading order.  The two
        diverge at finite strains (Y0/(2G) >= 0.05).
        """
        a: float = math.log(rho / self.rho_0)
        return (4.0 * self.G) / (3.0 * self.rho_0) * (1.0 - (1.0 + a) * math.exp(-a))

    # -- core algorithm ----------------------------------------------------

    def _solve_yield_and_precursor(self) -> None:
        """Compute yield state (rho_Y, e_Y, P_Y) and precursor (U_se, v_Y).

        Implements steps 2--5 of the solver algorithm.
        """
        # Step 2 -- yield density
        self.rho_Y = self.rho_0 * math.exp(self.Y_0 / (2.0 * self.G))
        assert self.rho_Y > self.rho_0, (
            f"Yield density rho_Y={self.rho_Y} must exceed rho_0={self.rho_0}"
        )

        # Step 3 -- yield energy (root-finding)
        e_el_Y: float = self._elastic_energy_work_integral(self.rho_Y)
        e_offset: float = e_el_Y if self.energy_split else 0.0

        e_lo: float = self.e_initial
        e_hi: float = self.e_initial + 10.0 * max(
            self.C_0**2, self.v_piston**2, self.Y_0 / self.rho_0,
        )

        def f_yield(e_var: float) -> float:
            P = self._eos(e_var, self.rho_Y)
            return (
                e_var + e_offset
                - self.e_initial
                - 0.5 / (self.rho_Y * self.rho_0)
                * (P + 2.0 / 3.0 * self.Y_0)
                * (self.rho_Y - self.rho_0)
            )

        e_var_Y: float = brentq(f_yield, e_lo, e_hi)

        if self.energy_split:
            self.e_th_Y = e_var_Y
            self.e_el_Y = e_el_Y
            self.e_Y = e_var_Y + e_el_Y
        else:
            self.e_Y = e_var_Y

        self.P_Y = self._eos(e_var_Y, self.rho_Y)

        assert math.isfinite(self.P_Y), f"Non-finite P_Y={self.P_Y}"
        assert math.isfinite(self.e_Y), f"Non-finite e_Y={self.e_Y}"

        # Step 4 -- elastic precursor speed
        self.U_se = math.sqrt(
            self.rho_Y * (self.P_Y + 2.0 / 3.0 * self.Y_0)
            / (self.rho_0 * (self.rho_Y - self.rho_0))
        )

        assert self.U_se > 0, f"U_se={self.U_se} must be positive"
        assert math.isfinite(self.U_se), f"Non-finite U_se={self.U_se}"

        # Step 5 -- elastic particle velocity
        self.v_Y = (self.rho_Y - self.rho_0) / self.rho_Y * self.U_se

        assert math.isfinite(self.v_Y), f"Non-finite v_Y={self.v_Y}"

        self.Sx_precursor = -2.0 / 3.0 * self.Y_0

    def _solve_elastic_wave(self) -> None:
        """Solve the elastic shock via 0-to-E Rankine-Hugoniot jump.

        When ``v_piston <= v_Y`` the material does not reach yield.  The
        elastic wave is a shock from state 0 to state E where:

            S_x^E = 4/3 G ln(rho_0 / rho_E)   (sub-yield deviatoric stress)
            rho_E = rho_0 U_se / (U_se - v_piston)
            P_E   = rho_0 U_se v_piston + S_x^E
            e_E   = e_0 + 1/(2 rho_0 rho_E)(P_E - S_x^E)(rho_E - rho_0)

        The elastic energy is the work integral of deviatoric power
        S:D^el, computed via ``_elastic_energy_work_integral(rho)``.

        After solving, the yield-point attributes (rho_Y, P_Y, e_Y, U_se,
        v_Y) are overwritten with the actual sub-yield elastic state.
        """
        U_se_yield = self.U_se

        def f_elastic(U_se: float) -> float:
            rho_E = self.rho_0 * U_se / (U_se - self.v_piston)
            Sx_E = 4.0 / 3.0 * self.G * math.log(self.rho_0 / rho_E)
            P_E = self.rho_0 * U_se * self.v_piston + Sx_E
            e_E = (
                self.e_initial
                + 0.5 / (self.rho_0 * rho_E)
                * (P_E - Sx_E) * (rho_E - self.rho_0)
            )
            if self.energy_split:
                e_eos = e_E - self._elastic_energy_work_integral(rho_E)
            else:
                e_eos = e_E
            return P_E - self._eos(e_eos, rho_E)

        roots = _scan_sign_change_roots(
            f_elastic,
            self.v_piston * 1.001,
            max(U_se_yield * 2.0, self.v_piston * 10.0),
            samples=600,
        )
        if not roots:
            raise RuntimeError(
                "No root found for the elastic shock residual.  "
                f"v_piston={self.v_piston}, U_se_yield={U_se_yield}."
            )
        self.U_se = max(roots)

        self.rho_Y = self.rho_0 * self.U_se / (self.U_se - self.v_piston)
        self.Sx_precursor = (
            4.0 / 3.0 * self.G * math.log(self.rho_0 / self.rho_Y)
        )
        self.P_Y = self.rho_0 * self.U_se * self.v_piston + self.Sx_precursor
        self.e_Y = (
            self.e_initial
            + 0.5 / (self.rho_0 * self.rho_Y)
            * (self.P_Y - self.Sx_precursor) * (self.rho_Y - self.rho_0)
        )
        self.v_Y = self.v_piston

        if self.energy_split:
            self.e_el_Y = self._elastic_energy_work_integral(self.rho_Y)
            self.e_th_Y = self.e_Y - self.e_el_Y

        assert math.isfinite(self.U_se), f"Non-finite U_se={self.U_se}"
        assert self.U_se > self.v_piston, (
            f"U_se={self.U_se} must exceed v_piston={self.v_piston}"
        )
        assert self.rho_Y > self.rho_0, (
            f"rho_E={self.rho_Y} must exceed rho_0={self.rho_0}"
        )
        assert math.isfinite(self.P_Y), f"Non-finite P_E={self.P_Y}"
        assert math.isfinite(self.e_Y), f"Non-finite e_E={self.e_Y}"

    def _classify_and_solve_shock(self) -> None:
        """Classify regime by comparing the direct-shock speed to U_se.

        First solves the 0-to-2 Rankine-Hugoniot jump to obtain U_s_02.
        If U_s_02 > U_se the plastic shock overruns the elastic precursor
        and a single plastic wave results.  Otherwise the two-wave structure
        is viable: the Y-to-2 residual is scanned below U_se and the
        largest root (weak-shock, nearest U_se) is selected.
        """
        # -- 0-to-2 residual ------------------------------------------------
        def f_shock_02(U_s: float) -> float:
            rho_2 = self.rho_0 * U_s / (U_s - self.v_piston)
            P_2 = (self.rho_0 * U_s * self.v_piston
                   - 2.0 / 3.0 * self.Y_0)
            e_2 = (
                self.e_initial
                + 0.5 / (self.rho_0 * rho_2)
                * (P_2 + 2.0 / 3.0 * self.Y_0) * (rho_2 - self.rho_0)
            )
            e_eos = (e_2 - self.e_el_Y) if self.energy_split else e_2
            return P_2 - self._eos(e_eos, rho_2)

        roots_02 = _scan_sign_change_roots(
            f_shock_02,
            self.v_piston * 1.001,
            max(2.0 * self.U_se, 5.0 * self.v_piston),
            samples=600,
        )
        if not roots_02:
            raise RuntimeError(
                "No shock root found in the 0-to-2 Rankine-Hugoniot "
                f"residual.  v_piston={self.v_piston}, "
                f"U_se={self.U_se}."
            )
        U_s_02 = max(roots_02)

        if U_s_02 > self.U_se:
            # -- Plastic wave: shock overruns precursor ----------------------
            self.wave_structure = WaveStructure.PLASTIC_WAVE
            self.U_s = U_s_02
            self._solve_shocked_state_02(self.U_s)
        else:
            # -- Two-wave: solve Y-to-2 below U_se ---------------------------
            self.wave_structure = WaveStructure.PLASTIC_SHOCK_AND_ELASTIC_WAVE

            def f_shock_y2(U_s: float) -> float:
                P_2 = (self.P_Y + self.rho_Y
                       * (U_s - self.v_Y) * (self.v_piston - self.v_Y))
                rho_2 = self.rho_Y * (U_s - self.v_Y) / (U_s - self.v_piston)
                e_2 = (
                    self.e_Y
                    + 0.5 / (self.rho_Y * rho_2)
                    * (self.P_Y + P_2 + 4.0 / 3.0 * self.Y_0)
                    * (rho_2 - self.rho_Y)
                )
                e_eos = (e_2 - self.e_el_Y) if self.energy_split else e_2
                return P_2 - self._eos(e_eos, rho_2)

            roots_y2 = _scan_sign_change_roots(
                f_shock_y2,
                self.v_piston * 1.001,
                self.U_se * (1.0 - 1.0e-10),
                samples=600,
            )
            if not roots_y2:
                raise RuntimeError(
                    "No Y-to-2 shock root found below U_se.  "
                    f"v_piston={self.v_piston}, U_se={self.U_se}."
                )
            self.U_s = max(roots_y2)

            assert self.U_s < self.U_se, (
                f"Two-wave: U_s={self.U_s} must be less than U_se={self.U_se}"
            )

            self._solve_shocked_state_y2(self.U_s)

    def _solve_shocked_state_y2(self, U_s: float) -> None:
        """Compute state 2 from Y-to-2 Rankine-Hugoniot jump.

        Both the Y and _2 states are at yield (S_x = -2/3 Y_0), so the
        deviatoric stress cancels in momentum and doubles in energy:

            rho_2 = rho_Y * (U_s - v_Y) / (U_s - v_piston)
            P_2   = P_Y + rho_Y * (U_s - v_Y) * (v_piston - v_Y)
            e_2   = e_Y + 1/(2*rho_Y*rho_2) * (P_Y + P_2 + 4/3*Y_0)
                        * (rho_2 - rho_Y)
        """
        self.rho_2 = (self.rho_Y * (U_s - self.v_Y)
                      / (U_s - self.v_piston))
        self.P_2 = (self.P_Y + self.rho_Y
                    * (U_s - self.v_Y) * (self.v_piston - self.v_Y))
        self.e_2 = (
            self.e_Y
            + 0.5 / (self.rho_Y * self.rho_2)
            * (self.P_Y + self.P_2 + 4.0 / 3.0 * self.Y_0)
            * (self.rho_2 - self.rho_Y)
        )

        if self.energy_split:
            self.e_el_2 = self.e_el_Y
            self.e_th_2 = self.e_2 - self.e_el_2

        assert math.isfinite(self.P_2), f"Non-finite P_2={self.P_2}"
        assert self.rho_2 > 0 and math.isfinite(self.rho_2), (
            f"Invalid rho_2={self.rho_2}"
        )
        assert math.isfinite(self.e_2), f"Non-finite e_2={self.e_2}"

    def _solve_shocked_state_02(self, U_s: float) -> None:
        """Compute state 2 from 0-to-2 Rankine-Hugoniot jump.

        State 0 has zero deviatoric stress (sigma_0 ~ 0), state 2 is at yield
        (S_x = -2/3 Y_0):

            rho_2 = rho_0 * U_s / (U_s - v_piston)
            P_2   = rho_0 * U_s * v_piston - 2/3 * Y_0
            e_2   = e_0 + 1/(2*rho_0*rho_2) * (P_2 + 2/3*Y_0)
                        * (rho_2 - rho_0)
        """
        self.rho_2 = self.rho_0 * U_s / (U_s - self.v_piston)
        self.P_2 = self.rho_0 * U_s * self.v_piston - 2.0 / 3.0 * self.Y_0
        self.e_2 = (
            self.e_initial
            + 0.5 / (self.rho_0 * self.rho_2)
            * (self.P_2 + 2.0 / 3.0 * self.Y_0) * (self.rho_2 - self.rho_0)
        )

        if self.energy_split:
            self.e_el_2 = self.e_el_Y
            self.e_th_2 = self.e_2 - self.e_el_2

        assert math.isfinite(self.P_2), f"Non-finite P_2={self.P_2}"
        assert self.rho_2 > 0 and math.isfinite(self.rho_2), f"Invalid rho_2={self.rho_2}"
        assert math.isfinite(self.e_2), f"Non-finite e_2={self.e_2}"

    def _solve_wave_structure(self) -> None:
        """Compute all intermediate and final wave-structure quantities.

        Delegates to :meth:`_solve_yield_and_precursor` for the elastic state,
        then either classifies the regime as elastic-only or invokes
        :meth:`_classify_and_solve_shock` for the plastic-shock solve.
        """
        self._solve_yield_and_precursor()

        if self.v_piston <= self.v_Y:
            self._solve_elastic_wave()
            self.wave_structure = WaveStructure.ELASTIC_WAVE
        else:
            self._classify_and_solve_shock()

        self._validate()

    def _validate(self) -> None:
        """Check physics invariants conditional on wave structure."""
        ws = self.wave_structure
        assert ws is not None

        if ws == WaveStructure.ELASTIC_WAVE:
            for name in ("rho_Y", "e_Y", "P_Y", "U_se", "v_Y"):
                val = getattr(self, name)
                if val is None or not math.isfinite(val):
                    raise ValueError(
                        f"Non-finite or missing quantity: {name}={val}"
                    )
            tol = 1.0e-6 * self.Y_0
            if abs(self.Sx_precursor) > 2.0 / 3.0 * self.Y_0 + tol:
                raise ValueError(
                    f"Elastic wave: |Sx_precursor|={abs(self.Sx_precursor):.6e}"
                    f" exceeds sub-yield limit 2/3*Y_0={2.0/3.0*self.Y_0:.6e}"
                )

        elif ws == WaveStructure.PLASTIC_SHOCK_AND_ELASTIC_WAVE:
            if not self.U_se > self.U_s:
                raise ValueError(
                    f"Two-wave: U_se={self.U_se} must exceed U_s={self.U_s}."
                )
            if not (self.rho_2 > self.rho_Y > self.rho_0):
                raise ValueError(
                    f"Two-wave: monotonic compression violated: "
                    f"rho_0={self.rho_0}, rho_Y={self.rho_Y}, "
                    f"rho_2={self.rho_2}."
                )
            for name in ("rho_Y", "e_Y", "P_Y", "U_se", "v_Y",
                         "U_s", "rho_2", "P_2", "e_2"):
                val = getattr(self, name)
                if val is None or not math.isfinite(val):
                    raise ValueError(
                        f"Non-finite or missing quantity: {name}={val}"
                    )

        elif ws == WaveStructure.PLASTIC_WAVE:
            if not self.rho_2 > self.rho_0:
                raise ValueError(
                    f"Plastic wave: rho_2={self.rho_2} must exceed "
                    f"rho_0={self.rho_0}."
                )
            for name in ("U_s", "rho_2", "P_2", "e_2"):
                val = getattr(self, name)
                if val is None or not math.isfinite(val):
                    raise ValueError(
                        f"Non-finite or missing quantity: {name}={val}"
                    )

    # -- public API --------------------------------------------------------

    def solve(
        self,
        t: float,
        x: ArrayLike,
    ) -> dict[str, object]:
        """Evaluate the piecewise-constant analytic solution.

        Parameters
        ----------
        t : float
            Time at which to evaluate the solution.  Must be non-negative.
        x : array-like of floats
            Spatial positions at which to evaluate the solution.  All values
            must be non-negative.

        Returns
        -------
        dict
            Keys and values:

            * ``"density"``  -- numpy array, density at each *x_i*.
            * ``"pressure"`` -- numpy array, pressure at each *x_i*.
            * ``"velocity"`` -- numpy array, velocity at each *x_i*.
            * ``"energy"``   -- numpy array, specific internal energy at
              each *x_i* (total energy, including elastic contribution).
            * ``"Sx"``       -- numpy array, deviatoric stress component
              *S_x* at each *x_i*.
            * ``"stress"``   -- numpy array, total axial stress
              *sigma_x = S_x - P* at each *x_i*.
            * ``"shock_location"``             -- float or ``np.nan``,
              plastic shock position *U_s * t*.  ``np.nan`` when no plastic
              shock exists (elastic-only regime).
            * ``"elastic_precursor_location"`` -- float or ``np.nan``,
              elastic precursor position *U_se * t*.  ``np.nan`` when the
              elastic precursor is absent (one-wave plastic regime).
            * ``"piston_location"`` -- float, piston face position
              *v_piston * t*.
            * ``"wave_structure"`` -- :class:`WaveStructure` enum value.

            When ``energy_split=True``, two additional keys are present:

            * ``"e_thermal"`` -- numpy array, thermal specific energy.
            * ``"e_elastic"`` -- numpy array, elastic specific energy
              (work-integral formulation; frozen at yield value for
              plastic states).

        Notes
        -----
        All field arrays are set to ``np.nan`` at positions behind the
        piston face (*x < v_piston * t*), since the solution is not
        defined there.

        Sign convention
        ---------------
        Compression is **negative**: *S_x < 0* and *sigma_x < 0* in
        compressed regions.  This follows the convention in the reference
        derivation.
        """
        if t < 0:
            raise ValueError(f"Time must be non-negative, got t={t}.")

        x_arr = np.asarray(x, dtype=np.float64)

        if np.any(x_arr < 0):
            raise ValueError("All spatial positions x must be non-negative.")

        ws = self.wave_structure

        density = np.full_like(x_arr, self.rho_0)
        pressure = np.full_like(x_arr, self.P_initial)
        velocity = np.zeros_like(x_arr)
        energy = np.full_like(x_arr, self.e_initial)
        Sx = np.zeros_like(x_arr)

        if ws == WaveStructure.ELASTIC_WAVE:
            shock_loc: float = np.nan
            elastic_loc: float = self.U_se * t

            elastic = x_arr <= elastic_loc
            density[elastic] = self.rho_Y
            pressure[elastic] = self.P_Y
            velocity[elastic] = self.v_Y
            energy[elastic] = self.e_Y
            Sx[elastic] = self.Sx_precursor

        elif ws == WaveStructure.PLASTIC_SHOCK_AND_ELASTIC_WAVE:
            shock_loc = self.U_s * t
            elastic_loc = self.U_se * t

            shocked = x_arr < shock_loc
            elastic = (x_arr >= shock_loc) & (x_arr <= elastic_loc)

            density[elastic] = self.rho_Y
            pressure[elastic] = self.P_Y
            velocity[elastic] = self.v_Y
            energy[elastic] = self.e_Y
            Sx[elastic] = self.Sx_precursor

            density[shocked] = self.rho_2
            pressure[shocked] = self.P_2
            velocity[shocked] = self.v_piston
            energy[shocked] = self.e_2
            Sx[shocked] = -2.0 / 3.0 * self.Y_0

        elif ws == WaveStructure.PLASTIC_WAVE:
            shock_loc = self.U_s * t
            elastic_loc = np.nan

            shocked = x_arr < shock_loc
            density[shocked] = self.rho_2
            pressure[shocked] = self.P_2
            velocity[shocked] = self.v_piston
            energy[shocked] = self.e_2
            Sx[shocked] = -2.0 / 3.0 * self.Y_0

        stress = Sx - pressure

        piston_loc: float = self.v_piston * t
        behind_piston = x_arr < piston_loc
        for arr in (density, pressure, velocity, energy, Sx, stress):
            arr[behind_piston] = np.nan

        result: dict[str, object] = {
            "density": density,
            "pressure": pressure,
            "velocity": velocity,
            "energy": energy,
            "Sx": Sx,
            "stress": stress,
            "shock_location": shock_loc,
            "elastic_precursor_location": elastic_loc,
            "piston_location": piston_loc,
            "wave_structure": ws,
        }

        if self.energy_split:
            e_thermal = np.full_like(x_arr, self.e_initial)
            e_elastic = np.zeros_like(x_arr)

            if ws == WaveStructure.ELASTIC_WAVE:
                e_thermal[elastic] = self.e_th_Y
                e_elastic[elastic] = self.e_el_Y

            elif ws == WaveStructure.PLASTIC_SHOCK_AND_ELASTIC_WAVE:
                e_thermal[elastic] = self.e_th_Y
                e_elastic[elastic] = self.e_el_Y
                e_thermal[shocked] = self.e_th_2
                e_elastic[shocked] = self.e_el_2

            elif ws == WaveStructure.PLASTIC_WAVE:
                e_thermal[shocked] = self.e_th_2
                e_elastic[shocked] = self.e_el_2

            e_thermal[behind_piston] = np.nan
            e_elastic[behind_piston] = np.nan
            result["e_thermal"] = e_thermal
            result["e_elastic"] = e_elastic

        return result
