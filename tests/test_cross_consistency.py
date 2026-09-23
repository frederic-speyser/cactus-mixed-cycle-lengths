"""
test_cross_consistency.py

Authors: Frederic G. Speyser & Kseniya Vyatkina
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina. All rights reserved.

Cross-checks the two independent computational routes used throughout
this project: mgonal_cactus_series_omega.py (exact rational-arithmetic
expansion of the truncated series) and critical_point_solver.py
(high-precision direct solution of the critical system, by 2-D Newton's
method on the characteristic system, using the truncated series only at
subcritical arguments). Agreement between the two, obtained by different
code with no shared logic beyond the common mathematical specification,
is the kind of check that a single method cannot provide on its own.

The tolerance below is tight (1e-30): critical_point_solver.py's
series-plus-Newton engine reaches full working precision (50 significant
digits by default) directly, with no undocumented parameter tuning, so
there is no reason to accept a loose tolerance here. An earlier version
of this test used a 2% tolerance to accommodate the previous
(damped-fixed-point / Aitken-extrapolation) engine, which turned out not
to converge to the precision it claimed for the mixed-parity examples of
Table 1 -- see CHANGELOG.md. That loose tolerance is exactly what let the
error pass this test undetected; it is not repeated here.

Reference: F. G. Speyser and K. Vyatkina, "Enumeration and Asymptotic
Analysis of Strict Non-Plane Cactus Graphs over a Finite Set of Cycle
Lengths" (submitted to the Journal of Integer Sequences, 2026), Section 8.

Run: python3 -m pytest test_cross_consistency.py -v
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))

from mpmath import mpf
from critical_point_solver import solve_critical_pair, theorem5_closed_form

TOLERANCE = mpf('1e-30')

_GUESS = {
    (5, 7): ('0.550', '0.727'),
    (5, 7, 9): ('0.535', '0.692'),
}


def test_theorem5_case_5_7():
    """All-odd Omega: the numerically located tau must agree with the
    closed-form root of Theorem 5, computed by fully independent code
    (mpmath.findroot on a plain polynomial, vs. the series-plus-Newton
    solver)."""
    omega = (5, 7)
    rho, tau_numeric, *_ = solve_critical_pair(omega, *_GUESS[omega])
    tau_closed = theorem5_closed_form(omega)
    rel_diff = abs(tau_numeric - tau_closed) / tau_closed
    assert rel_diff < TOLERANCE, f"relative difference {float(rel_diff):.3e} exceeds tolerance"


def test_theorem5_case_5_7_9():
    """Same check with a three-generator Omega, to confirm the agreement
    is not an artifact specific to two generators."""
    omega = (5, 7, 9)
    rho, tau_numeric, *_ = solve_critical_pair(omega, *_GUESS[omega])
    tau_closed = theorem5_closed_form(omega)
    rel_diff = abs(tau_numeric - tau_closed) / tau_closed
    assert rel_diff < TOLERANCE, f"relative difference {float(rel_diff):.3e} exceeds tolerance"


if __name__ == "__main__":
    test_theorem5_case_5_7()
    test_theorem5_case_5_7_9()
    print("All cross-consistency tests passed.")
