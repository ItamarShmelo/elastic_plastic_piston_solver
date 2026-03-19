"""Analytic solver for the 1D elastoplastic piston test problem.

Solves the problem of a piston driving into an elastoplastic medium with a
Mie-Gruneisen equation of state.  The solution consists of three piecewise-
constant regions separated by an elastic precursor wave and a plastic shock.

Reference
---------
See the derivation in ``README.md`` for the full notation and equations.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import brentq


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
# Solver class
# ---------------------------------------------------------------------------

class ElastoplasticPistonSolver:
    """Analytic solver for the 1D elastoplastic piston problem.

    The solver computes the three-region piecewise-constant solution produced
    when a piston drives into an elastoplastic medium at constant velocity.

    Parameters
    ----------
    rho_0 : float
        Reference density.
    C_0 : float
        Bulk sound speed at reference state.
    s : float
        Hugoniot slope coefficient (linear Us-Up relation).
    Gamma_0 : float
        Gruneisen parameter (assumed constant).
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
        energy is ``S:S / (4 rho G)``.  Default is ``False`` (original
        single-energy formulation).
    """

    def __init__(
        self,
        rho_0: float,
        C_0: float,
        s: float,
        Gamma_0: float,
        G: float,
        Y_0: float,
        e_initial: float,
        v_piston: float,
        energy_split: bool = False,
    ) -> None:
        # Store EOS parameters
        self.rho_0: float = rho_0
        self.C_0: float = C_0
        self.s: float = s
        self.Gamma_0: float = Gamma_0

        # Elastoplastic parameters
        self.G: float = G
        self.Y_0: float = Y_0

        # Initial / boundary conditions
        self.e_initial: float = e_initial
        self.v_piston: float = v_piston

        # Mode flag
        self.energy_split: bool = energy_split

        # Solved state variables (populated by _solve_wave_structure)
        self.rho_Y: float | None = None
        self.e_Y: float | None = None
        self.P_Y: float | None = None
        self.U_se: float | None = None
        self.v_Y: float | None = None
        self.U_s: float | None = None
        self.rho_2: float | None = None
        self.P_2: float | None = None
        self.e_2: float | None = None

        # Energy-split mode quantities (None when energy_split=False)
        self.e_th_Y: float | None = None
        self.e_el_Y: float | None = None
        self.e_th_2: float | None = None
        self.e_el_2: float | None = None

        # Initial-region pressure (full EOS at reference state)
        self.P_initial: float = self._eos(e_initial, rho_0)

        self._solve_wave_structure()

    # -- convenience wrappers ----------------------------------------------

    def _eos(self, e: float, rho: float) -> float:
        """Evaluate the Mie-Gruneisen EOS with stored parameters."""
        return _mie_gruneisen_pressure(
            e, rho, self.rho_0, self.C_0, self.s, self.Gamma_0,
        )

    # -- core algorithm ----------------------------------------------------

    def _solve_wave_structure(self) -> None:
        """Compute all intermediate and final wave-structure quantities.

        Implements steps 2-7 of the solver algorithm from the reference,
        including root-finding for e^Y and U_s.
        """
        rho_0 = self.rho_0
        C_0 = self.C_0
        Y_0 = self.Y_0
        G = self.G
        e_0 = self.e_initial
        v_piston = self.v_piston

        # Step 2 -- yield density
        rho_Y: float = rho_0 * math.exp(Y_0 / (2.0 * G))
        self.rho_Y = rho_Y

        # Step 3 -- yield energy (root-finding)
        two_thirds_Y0: float = 2.0 / 3.0 * Y_0
        four_thirds_Y0: float = 4.0 / 3.0 * Y_0
        e_el_Y: float = Y_0**2 / (6.0 * rho_Y * G)

        e_lo: float = e_0
        e_hi: float = e_0 + 10.0 * C_0**2

        if self.energy_split:
            def f_yield(e_th: float) -> float:
                P = self._eos(e_th, rho_Y)
                return (
                    e_th + e_el_Y
                    - e_0
                    - 0.5 / (rho_Y * rho_0) * (P + two_thirds_Y0) * (rho_Y - rho_0)
                )

            e_th_Y: float = brentq(f_yield, e_lo, e_hi)
            self.e_th_Y = e_th_Y
            self.e_el_Y = e_el_Y
            self.e_Y = e_th_Y + e_el_Y
            P_Y: float = self._eos(e_th_Y, rho_Y)
        else:
            def f_yield(e: float) -> float:
                P = self._eos(e, rho_Y)
                return (
                    e
                    - e_0
                    - 0.5 / (rho_Y * rho_0) * (P + two_thirds_Y0) * (rho_Y - rho_0)
                )

            e_Y: float = brentq(f_yield, e_lo, e_hi)
            self.e_Y = e_Y
            P_Y: float = self._eos(e_Y, rho_Y)

        self.P_Y = P_Y

        # Step 4 -- elastic precursor speed
        #   U_se^2 = rho_Y * (P_Y + 2/3 Y_0) / (rho_0 * (rho_Y - rho_0))
        U_se_sq: float = rho_Y * (P_Y + two_thirds_Y0) / (rho_0 * (rho_Y - rho_0))
        U_se: float = math.sqrt(U_se_sq)
        self.U_se = U_se

        # Step 5 -- elastic particle velocity
        v_Y: float = (rho_Y - rho_0) / rho_Y * U_se
        self.v_Y = v_Y

        # Check that the piston is fast enough to produce a plastic shock
        if v_piston <= v_Y:
            raise ValueError(
                f"Piston velocity v_piston={v_piston} must exceed the "
                f"elastic-limit particle velocity v_Y={v_Y:.4f} for a "
                f"plastic shock to form."
            )

        # Step 6 -- plastic shock speed (root-finding)
        #
        # The equation f_shock(U_s) = P_2(U_s) - P_eos(rho_2, e_2) can have
        # multiple roots.  The physical root is the one nearest U_se (moderate
        # compression), not the spurious low-velocity root (extreme
        # compression).  We search downward from U_se to find the bracket.
        if self.energy_split:
            e_th_Y_local = self.e_th_Y
            assert e_th_Y_local is not None

            def f_shock(U_s: float) -> float:
                v_2 = v_piston
                P_2 = P_Y + rho_Y * (U_s - v_Y) * (v_2 - v_Y)
                rho_2 = rho_Y * (U_s - v_Y) / (U_s - v_2)
                e_th_2 = (
                    e_th_Y_local
                    + e_el_Y
                    + 0.5 / (rho_Y * rho_2)
                    * (P_Y + P_2 + four_thirds_Y0)
                    * (rho_2 - rho_Y)
                    - Y_0**2 / (6.0 * rho_2 * G)
                )
                return P_2 - self._eos(e_th_2, rho_2)
        else:
            e_Y_local = self.e_Y
            assert e_Y_local is not None

            def f_shock(U_s: float) -> float:
                v_2 = v_piston
                P_2 = P_Y + rho_Y * (U_s - v_Y) * (v_2 - v_Y)
                rho_2 = rho_Y * (U_s - v_Y) / (U_s - v_2)
                e_2 = (
                    e_Y_local
                    + 0.5 / (rho_Y * rho_2)
                    * (P_Y + P_2 + four_thirds_Y0)
                    * (rho_2 - rho_Y)
                )
                return P_2 - self._eos(e_2, rho_2)

        # Bracket: search downward from U_se in geometric steps until
        # f_shock changes sign.  This finds the physical root (nearest U_se).
        U_s_hi: float = U_se * (1.0 - 1.0e-10)
        f_hi: float = f_shock(U_s_hi)
        U_s_lo: float = U_s_hi
        ratio: float = 0.9  # geometric step factor
        for _ in range(200):
            U_s_lo = max(U_s_lo * ratio, v_piston * 1.001)
            f_lo: float = f_shock(U_s_lo)
            if f_lo * f_hi < 0.0:
                break
            if U_s_lo <= v_piston * 1.001:
                raise RuntimeError(
                    "Could not bracket the plastic shock speed.  "
                    "Check input parameters."
                )
        else:
            raise RuntimeError(
                "Could not bracket the plastic shock speed after 200 "
                "iterations.  Check input parameters."
            )

        U_s: float = brentq(f_shock, U_s_lo, U_s_hi)
        self.U_s = U_s

        # Step 7 -- shocked-state quantities
        v_2 = v_piston
        P_2: float = P_Y + rho_Y * (U_s - v_Y) * (v_2 - v_Y)
        rho_2: float = rho_Y * (U_s - v_Y) / (U_s - v_2)

        if self.energy_split:
            assert self.e_th_Y is not None
            e_el_2: float = Y_0**2 / (6.0 * rho_2 * G)
            e_th_2: float = (
                self.e_th_Y
                + e_el_Y
                + 0.5 / (rho_Y * rho_2)
                * (P_Y + P_2 + four_thirds_Y0)
                * (rho_2 - rho_Y)
                - e_el_2
            )
            self.e_th_2 = e_th_2
            self.e_el_2 = e_el_2
            self.e_2 = e_th_2 + e_el_2
        else:
            self.e_2 = (
                self.e_Y
                + 0.5 / (rho_Y * rho_2)
                * (P_Y + P_2 + four_thirds_Y0)
                * (rho_2 - rho_Y)
            )

        self.P_2 = P_2
        self.rho_2 = rho_2

        # Step 8 -- physics invariant checks
        self._validate()

    def _validate(self) -> None:
        """Check physics invariants; raise ``ValueError`` on failure."""
        assert self.U_se is not None and self.U_s is not None
        assert self.rho_Y is not None and self.rho_2 is not None
        if not self.U_se > self.U_s:
            raise ValueError(
                f"Elastic precursor speed U_se={self.U_se} must exceed "
                f"plastic shock speed U_s={self.U_s}."
            )
        if not (self.rho_2 > self.rho_Y > self.rho_0):
            raise ValueError(
                f"Monotonic compression violated: "
                f"rho_0={self.rho_0}, rho_Y={self.rho_Y}, rho_2={self.rho_2}."
            )
        for name in ("rho_Y", "e_Y", "P_Y", "U_se", "v_Y", "U_s", "rho_2", "P_2", "e_2"):
            val = getattr(self, name)
            if not math.isfinite(val):
                raise ValueError(f"Non-finite solved quantity: {name}={val}")

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
            * ``"shock_location"``             -- float, plastic shock
              position *U_s * t*.
            * ``"elastic_precursor_location"`` -- float, elastic precursor
              position *U_se * t*.

            When ``energy_split=True``, two additional keys are present:

            * ``"e_thermal"`` -- numpy array, thermal specific energy.
            * ``"e_elastic"`` -- numpy array, elastic specific energy
              (*S:S / (4 rho G)*).

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

        assert self.U_s is not None and self.U_se is not None
        assert self.rho_Y is not None and self.e_Y is not None and self.P_Y is not None
        assert self.v_Y is not None
        assert self.rho_2 is not None and self.P_2 is not None and self.e_2 is not None

        shock_loc: float = self.U_s * t
        elastic_loc: float = self.U_se * t

        # Region masks
        shocked = x_arr < shock_loc
        elastic = (x_arr >= shock_loc) & (x_arr <= elastic_loc)
        # initial = x_arr > elastic_loc  (default)

        two_thirds_Y0: float = 2.0 / 3.0 * self.Y_0

        # Allocate output arrays initialised to the *initial* region values
        density = np.full_like(x_arr, self.rho_0)
        pressure = np.full_like(x_arr, self.P_initial)
        velocity = np.zeros_like(x_arr)
        energy = np.full_like(x_arr, self.e_initial)
        Sx = np.zeros_like(x_arr)

        # Elastic region
        density[elastic] = self.rho_Y
        pressure[elastic] = self.P_Y
        velocity[elastic] = self.v_Y
        energy[elastic] = self.e_Y
        Sx[elastic] = -two_thirds_Y0

        # Shocked region
        density[shocked] = self.rho_2
        pressure[shocked] = self.P_2
        velocity[shocked] = self.v_piston
        energy[shocked] = self.e_2
        Sx[shocked] = -two_thirds_Y0

        # Total axial stress: sigma_x = S_x - P
        stress = Sx - pressure

        result: dict[str, object] = {
            "density": density,
            "pressure": pressure,
            "velocity": velocity,
            "energy": energy,
            "Sx": Sx,
            "stress": stress,
            "shock_location": shock_loc,
            "elastic_precursor_location": elastic_loc,
        }

        if self.energy_split:
            assert self.e_th_Y is not None and self.e_el_Y is not None
            assert self.e_th_2 is not None and self.e_el_2 is not None
            e_thermal = np.full_like(x_arr, self.e_initial)
            e_elastic = np.zeros_like(x_arr)

            e_thermal[elastic] = self.e_th_Y
            e_elastic[elastic] = self.e_el_Y

            e_thermal[shocked] = self.e_th_2
            e_elastic[shocked] = self.e_el_2

            result["e_thermal"] = e_thermal
            result["e_elastic"] = e_elastic

        return result
