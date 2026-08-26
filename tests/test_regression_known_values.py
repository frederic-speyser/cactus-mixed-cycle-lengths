"""
test_regression_known_values.py

Regression tests anchoring critical_point_solver.py against numerical
values already published in the companion paper [1] -- not against
values computed by this project's own code, so that a bug shared by
both branches of K_C_term (odd and even m) cannot silently pass
undetected by agreeing with itself.

  - test_omega5_matches_ejc_theorem2: Omega = {5} is a singleton, odd,
    covered by the *odd* branch of K_C_term; tau_5 = 2^(-1/4) is an
    exact closed form (Theorem 2 of [1]), used at full solver precision.
  - test_omega6_matches_ejc_table3: Omega = {6} is a singleton, even,
    covered by the *even* branch of K_C_term (the one exercised by
    Proposition 6 for mixed Omega, and otherwise untested against any
    external ground truth in this project). (rho_6, tau_6) are the
    numerical values reported in Table 3 of [1], to 6 significant
    digits; the tolerance below is set accordingly, looser than the
    exact odd-case check.

Reference:
  [1] F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
      Non-Plane m-Gonal Cactus Graphs via Split-Decomposition", Theorem 2
      (tau_5) and Table 3 (rho_6, tau_6).


Author: Frederic G. Speyser
Run: python3 -m pytest test_regression_known_values.py -v
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))

from mpmath import mpf
from critical_point_solver import solve_rho, solve_tau, theorem5_closed_form

TAU_5_EJC_THEOREM2 = mpf(2) ** (mpf(-1) / 4)   # exact, = 2^(-1/4)
RHO_6_EJC_TABLE3 = mpf('0.633235')             # 6 significant digits, [1] Table 3
TAU_6_EJC_TABLE3 = mpf('0.821008')             # 6 significant digits, [1] Table 3


def test_omega5_matches_ejc_theorem2():
    """Odd branch of K_C_term: Theorem 5's closed form, specialized to a
    singleton, must equal the exact value 2^(-1/4) of [1, Theorem 2]."""
    tau = theorem5_closed_form((5,))
    assert abs(tau - TAU_5_EJC_THEOREM2) < mpf('1e-15')


def test_omega5_solver_matches_ejc_theorem2():
    """The general-purpose numerical solver (not just the closed-form
    shortcut) must also recover tau_5, within the precision expected of
    its default parameters."""
    rho = solve_rho((5,), mpf('0.6'), mpf('0.61'))
    tau, _ = solve_tau((5,), rho)
    assert abs(tau - TAU_5_EJC_THEOREM2) / TAU_5_EJC_THEOREM2 < mpf('0.01')


def test_omega6_matches_ejc_table3():
    """Even branch of K_C_term (the branch otherwise exercised only by
    Proposition 6 for mixed Omega, with no other external anchor in this
    project): the solver must reproduce the published numerical values
    for the singleton even case Omega = {6}."""
    rho = solve_rho((6,), mpf('0.62'), mpf('0.65'))
    tau, _ = solve_tau((6,), rho)
    assert abs(rho - RHO_6_EJC_TABLE3) / RHO_6_EJC_TABLE3 < mpf('0.01')
    assert abs(tau - TAU_6_EJC_TABLE3) / TAU_6_EJC_TABLE3 < mpf('0.01')


if __name__ == "__main__":
    test_omega5_matches_ejc_theorem2()
    test_omega5_solver_matches_ejc_theorem2()
    test_omega6_matches_ejc_table3()
    print("All regression tests against published values passed.")
