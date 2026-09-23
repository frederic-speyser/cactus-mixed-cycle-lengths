"""
test_regression_known_values.py

Authors: Frederic G. Speyser & Kseniya Vyatkina
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina. All rights reserved.

Regression tests anchoring critical_point_solver.py against numerical
values already published in the companion paper [1] -- not against
values computed by this project's own code, so that a bug shared by both
branches of the kernel (odd and even m) cannot silently pass undetected
by agreeing with itself.

  - test_omega5_matches_ejc_theorem2: Omega = {5} is a singleton, odd,
    covered by the *odd* branch of the kernel; tau_5 = 2^(-1/4) is an
    exact closed form (Theorem 2 of [1]), used at full solver precision.
  - test_omega6_matches_ejc_table3: Omega = {6} is a singleton, even,
    covered by the *even* branch of the kernel (the one exercised by
    Proposition 6 for mixed Omega, and otherwise untested against any
    external ground truth in this project). (rho_6, tau_6) are the
    numerical values reported in Table 3 of [1], to 6 significant
    digits, so the tolerance below (1e-6 relative) is set to that
    precision -- not loosened further, unlike the previous version of
    this test (see CHANGELOG.md).

Reference:
  [1] F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
      Non-Plane m-Gonal Cactus Graphs via Split-Decomposition", Theorem 2
      (tau_5) and Table 3 (rho_6, tau_6).

Run: python3 -m pytest test_regression_known_values.py -v
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))

from mpmath import mpf
from critical_point_solver import solve_critical_pair, theorem5_closed_form

TAU_5_EJC_THEOREM2 = mpf(2) ** (mpf(-1) / 4)   # exact, = 2^(-1/4)
RHO_6_EJC_TABLE3 = mpf('0.633235')             # 6 significant digits, [1] Table 3
TAU_6_EJC_TABLE3 = mpf('0.821008')             # 6 significant digits, [1] Table 3

EJC_TOLERANCE = mpf('1e-6')  # matches the 6-significant-digit precision of [1]'s table


def test_omega5_matches_ejc_theorem2():
    """Odd branch of the kernel: Theorem 5's closed form, specialized to
    a singleton, must equal the exact value 2^(-1/4) of [1, Theorem 2]."""
    tau = theorem5_closed_form((5,))
    assert abs(tau - TAU_5_EJC_THEOREM2) < mpf('1e-15')


def test_omega5_solver_matches_ejc_theorem2():
    """The general-purpose numerical solver (not just the closed-form
    shortcut) must also recover tau_5, at full solver precision."""
    rho, tau, *_ = solve_critical_pair((5,), '0.604', '0.840')
    assert abs(tau - TAU_5_EJC_THEOREM2) / TAU_5_EJC_THEOREM2 < mpf('1e-30')


def test_omega6_matches_ejc_table3():
    """Even branch of the kernel (the branch otherwise exercised only by
    Proposition 6 for mixed Omega, with no other external anchor in this
    project): the solver must reproduce the published numerical values
    for the singleton even case Omega = {6}, to the precision [1]
    reports them."""
    rho, tau, *_ = solve_critical_pair((6,), '0.633', '0.821')
    assert abs(rho - RHO_6_EJC_TABLE3) / RHO_6_EJC_TABLE3 < EJC_TOLERANCE
    assert abs(tau - TAU_6_EJC_TABLE3) / TAU_6_EJC_TABLE3 < EJC_TOLERANCE


if __name__ == "__main__":
    test_omega5_matches_ejc_theorem2()
    test_omega5_solver_matches_ejc_theorem2()
    test_omega6_matches_ejc_table3()
    print("All regression tests against published values passed.")
