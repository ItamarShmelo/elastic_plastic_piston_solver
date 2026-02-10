"""Analytic solver for the 1D elastoplastic piston test problem.

Solves the problem of a piston driving into an elastoplastic medium with a
Mie-Gruneisen equation of state.  The solution consists of three piecewise-
constant regions separated by an elastic precursor wave and a plastic shock.

Reference
---------
See ``references/elastoplastic_piston_test_problem.tex`` for the full
derivation and notation.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import numpy.typing as npt
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
    P_0 : float
        Reference pressure offset.
    G : float
        Shear modulus.
    Y_0 : float
        Yield strength.
    e_initial : float
        Initial specific internal energy.
    v_piston : float
        Piston velocity (positive, in the +x direction).
    """

    def __init__(
        self,
        rho_0: float,
        C_0: float,
        s: float,
        Gamma_0: float,
        P_0: float,
        G: float,
        Y_0: float,
        e_initial: float,
        v_piston: float,
    ) -> None:
        # Store EOS parameters
        self.rho_0: float = rho_0
        self.C_0: float = C_0
        self.s: float = s
        self.Gamma_0: float = Gamma_0
        self.P_0: float = P_0

        # Elastoplastic parameters
        self.G: float = G
        self.Y_0: float = Y_0

        # Initial / boundary conditions
        self.e_initial: float = e_initial
        self.v_piston: float = v_piston

        # Solved state variables (populated by _solve_wave_structure)
        self.rho_Y: float = 0.0
        self.e_Y: float = 0.0
        self.P_Y: float = 0.0
        self.U_se: float = 0.0
        self.v_Y: float = 0.0
        self.U_s: float = 0.0
        self.rho_2: float = 0.0
        self.P_2: float = 0.0
        self.e_2: float = 0.0

        # Initial-region pressure (full EOS at reference state)
        self.P_initial: float = self._eos(e_initial, rho_0)

        self._solve_wave_structure()

    # -- convenience wrappers ----------------------------------------------

    def _eos(self, e: float, rho: float) -> float:
        """Evaluate the Mie-Gruneisen EOS with stored parameters."""
        return _mie_gruneisen_pressure(
            e, rho, self.rho_0, self.C_0, self.s, self.Gamma_0, self.P_0,
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

        def f_yield(e: float) -> float:
            P = self._eos(e, rho_Y)
            return (
                e
                - e_0
                - 0.5 / (rho_Y * rho_0) * (P + two_thirds_Y0) * (rho_Y - rho_0)
            )

        e_lo: float = e_0
        e_hi: float = e_0 + 10.0 * C_0**2
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
        def f_shock(U_s: float) -> float:
            v_2 = v_piston
            P_2 = P_Y + rho_Y * (U_s - v_Y) * (v_2 - v_Y)
            rho_2 = rho_Y * (U_s - v_Y) / (U_s - v_2)
            e_2 = (
                e_Y
                + 0.5 / (rho_Y * rho_2)
                * (P_Y + P_2 + four_thirds_Y0)
                * (rho_2 - rho_Y)
            )
            P_eos = self._eos(e_2, rho_2)
            return P_2 - P_eos

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
        e_2: float = (
            e_Y
            + 0.5 / (rho_Y * rho_2)
            * (P_Y + P_2 + four_thirds_Y0)
            * (rho_2 - rho_Y)
        )
        self.P_2 = P_2
        self.rho_2 = rho_2
        self.e_2 = e_2

        # Step 8 -- physics invariant checks
        self._validate()

    def _validate(self) -> None:
        """Check physics invariants; raise ``ValueError`` on failure."""
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
        x: npt.NDArray[np.floating[Any]],
    ) -> dict[str, object]:
        """Evaluate the piecewise-constant analytic solution.

        Parameters
        ----------
        t : float
            Time at which to evaluate the solution.
        x : numpy array of floats
            Spatial positions at which to evaluate the solution.

        Returns
        -------
        dict
            Keys and values:

            * ``"density"``  -- numpy array, density at each *x_i*.
            * ``"pressure"`` -- numpy array, pressure at each *x_i*.
            * ``"velocity"`` -- numpy array, velocity at each *x_i*.
            * ``"energy"``   -- numpy array, specific internal energy at
              each *x_i*.
            * ``"Sx"``       -- numpy array, deviatoric stress component
              *S_x* at each *x_i*.
            * ``"stress"``   -- numpy array, total axial stress
              *sigma_x = S_x - P* at each *x_i*.
            * ``"shock_location"``             -- float, plastic shock
              position *U_s * t*.
            * ``"elastic_precursor_location"`` -- float, elastic precursor
              position *U_se * t*.

        Sign convention
        ---------------
        Compression is **negative**: *S_x < 0* and *sigma_x < 0* in
        compressed regions.  This follows the convention in the reference
        derivation.
        """
        shock_loc: float = self.U_s * t
        elastic_loc: float = self.U_se * t

        x_arr: npt.NDArray[np.floating[Any]] = np.asarray(x, dtype=np.float64)

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

        return {
            "density": density,
            "pressure": pressure,
            "velocity": velocity,
            "energy": energy,
            "Sx": Sx,
            "stress": stress,
            "shock_location": shock_loc,
            "elastic_precursor_location": elastic_loc,
        }
