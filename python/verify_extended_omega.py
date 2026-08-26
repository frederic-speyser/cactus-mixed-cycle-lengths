"""
verify_extended_omega.py

Extended independent verification on a batch of finite sets Omega that do
NOT appear anywhere in Table 1 of [2] (the four "representative examples"
printed in the paper: {5,6}, {5,7}, {5,7,9}, {5,6,7}). The point of this
script is to check that the theorems and the numerical machinery of [2]
are not somehow implicitly tuned to those four cases: everything checked
here is proved in [2] for arbitrary finite Omega, so it should hold on
any Omega we throw at it, and this script throws eight more at it, chosen
to cover the same qualitative combinations as Table 1 (both sizes
|Omega| = 2 and |Omega| = 3, both the all-odd and the mixed-parity
regimes), plus two cases neither Table 1 nor the original six-Omega batch
covered:

    {7,11}       odd,   |Omega|=2   (untabulated)
    {6,7}        mixed, |Omega|=2   (untabulated)
    {5,8}        mixed, |Omega|=2   (untabulated)
    {7,9,11}     odd,   |Omega|=3   (untabulated)
    {5,9,13}     odd,   |Omega|=3   (untabulated; also used by
                                      test_theorem_a_untabulated_omega.py
                                      for Theorem 3 alone)
    {5,6,8}      mixed, |Omega|=3   (untabulated; one even length)
    {5,7,9,11}   odd,   |Omega|=4   (untabulated; Table 1 and the six
                                      Omega above never go past |Omega|=3)
    {6,8,9}      mixed, |Omega|=3   (untabulated; TWO even lengths at
                                      once -- every mixed case above,
                                      and every mixed case in Table 1,
                                      has exactly one even length, so
                                      Proposition 6's sum over Omega_even
                                      was never exercised with more than
                                      one term before this case)

For each Omega in this list, three independent things are checked:

  1. Theorem 3 (exact support characterization): the zero/nonzero
     pattern of the exact series from mgonal_cactus_series_omega.py is
     compared against a direct breadth-first search over the numerical
     semigroup Gamma_Omega, exactly as in
     test_theorem_a_untabulated_omega.py (reproduced here so this script
     is self-contained).

  2. Two independent series implementations agree term by term: the
     main solver, mgonal_cactus_series_omega.py (Python,
     fractions.Fraction, hand-written convolutions), against
     verify_dissymmetry_omega.py (Python, sympy.Rational, a separately
     written set of convolutions) -- both the rooted series s(x) and the
     unrooted series G(x), for as many terms as both scripts compute by
     default.

  3. The critical pair located numerically by critical_point_solver.py
     is self-consistent with the relevant closed-form/identity of [2]:
     Theorem 5's closed form for the all-odd cases, Proposition 6's
     implicit equation for the mixed-parity cases -- exactly the
     consistency checks already used for the four tabulated examples in
     reproduce_table1.py, just run here on untabulated Omega instead.

This script does not by itself constitute a fourth, language-independent
check -- that role is filled by verify_pari_omega_extended.gp, which
runs the same batch of six Omega through PARI/GP's native power series
arithmetic and reports the same rooted/unrooted terms for direct
comparison against the output printed here.

Reference: F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths" (in
preparation for submission to the Journal of Integer Sequences, 2026).

Auteur : Frederic G. Speyser
Run: python3 verify_extended_omega.py
"""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

import mgonal_cactus_series_omega as om
import verify_dissymmetry_omega as vd
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(__file__))
import critical_point_solver as cps
from mpmath import mpf, nstr

UNTABULATED_OMEGA = [
    (7, 11),
    (6, 7),
    (5, 8),
    (7, 9, 11),
    (5, 9, 13),
    (5, 6, 8),
    (5, 7, 9, 11),
    (6, 8, 9),
]

BOUND = 40          # for the Theorem 3 support check
CROSS_N = 30         # truncation order for the two-implementation cross-check
CROSS_TERMS = 8      # number of leading non-zero terms compared

# Numeric critical-point settings. The damped fixed-point iteration converges
# more slowly for some of the Omega above (empirically, ones involving a
# larger cycle length such as 11 or 13) than for the four Omega tabulated in
# the paper, for which the defaults in critical_point_solver.py itself were
# adequate. The values below were chosen empirically (see CHANGELOG.md) so
# that the consistency check below lands within ACCEPT_TOL for every Omega
# in UNTABULATED_OMEGA in practical runtime. This is still well short of the
# much higher-precision, longer-running settings used to produce the
# seven-decimal-digit values printed in Table 1 of the paper (see the module
# docstring of critical_point_solver.py); it is enough to confirm there is
# no structural inconsistency, which is this script's only purpose here.
RHO_ITERS_PROBE = 1600
RHO_K = 42
RHO_ROUNDS = 24
TAU_ITERS = 7000
TAU_K = 55
ACCEPT_TOL = mpf('0.005')  # relative (Theorem 5) / absolute out of 2 (Prop. 6)


def semigroup_support_bfs(omega, bound):
    gens = [m - 1 for m in omega]
    reachable = {0}
    frontier = {0}
    while frontier:
        new_frontier = set()
        for t in frontier:
            for g in gens:
                nt = t + g
                if nt <= bound - 1 and nt not in reachable:
                    reachable.add(nt)
                    new_frontier.add(nt)
        frontier = new_frontier
    return {t + 1 for t in reachable}


def check_theorem3(omega):
    om.N = BOUND
    s, _ = om.solve_G(omega)
    actual_support = {n for n in range(1, BOUND + 1) if s[n] != 0}
    predicted_support = semigroup_support_bfs(omega, BOUND)
    ok = actual_support == predicted_support
    return ok, actual_support, predicted_support


def leading_terms_fraction(coeffs, count):
    out = []
    for n, c in enumerate(coeffs):
        if c == 0:
            continue
        assert c.denominator == 1, f"non-integer coefficient at x^{n}: {c}"
        out.append((n, int(c)))
        if len(out) == count:
            break
    return out


def check_cross_implementation(omega):
    om.N = CROSS_N
    s_main, G_main = om.solve_G(omega)
    s_main_terms = leading_terms_fraction(s_main, CROSS_TERMS)
    G_main_terms = leading_terms_fraction(G_main, CROSS_TERMS)

    s_vd, G_vd = vd.solve_G(omega, CROSS_N)
    s_vd_terms = vd.nonzero_terms(s_vd, CROSS_TERMS)
    G_vd_terms = vd.nonzero_terms(G_vd, CROSS_TERMS)

    ok = (s_main_terms == s_vd_terms) and (G_main_terms == G_vd_terms)
    return ok, s_main_terms, G_main_terms, s_vd_terms, G_vd_terms


def check_critical_pair(omega):
    all_odd = all(m % 2 == 1 for m in omega)
    lo, hi = mpf('0.30'), mpf('0.70')
    rho = cps.solve_rho(omega, lo, hi, rounds=RHO_ROUNDS,
                         iters_probe=RHO_ITERS_PROBE, K=RHO_K)
    tau, V = cps.solve_tau(omega, rho, iters=TAU_ITERS, K=TAU_K)
    if all_odd:
        tau_exact = cps.theorem5_closed_form(omega)
        rel_diff = abs(tau - tau_exact) / tau_exact
        ok = rel_diff < ACCEPT_TOL
        return ok, ('theorem5', rho, tau, tau_exact, rel_diff)
    else:
        lhs = cps.proposition6_check(omega, rho, tau, V)
        diff = abs(lhs - 2)
        ok = diff < ACCEPT_TOL
        return ok, ('prop6', rho, tau, lhs, diff)


def main():
    all_ok = True
    omega_list = UNTABULATED_OMEGA
    if len(sys.argv) > 1:
        # optional "start:end" slice, e.g. "0:3", to split a full run across
        # several invocations (each Omega's numeric critical-pair check can
        # take a minute or more at the precision settings above).
        start, end = (int(v) for v in sys.argv[1].split(':'))
        omega_list = UNTABULATED_OMEGA[start:end]
    print("=" * 78)
    print("EXTENDED VERIFICATION ON UNTABULATED OMEGA")
    print("(none of these appear in Table 1 of the paper)")
    print("=" * 78)
    for omega in omega_list:
        print()
        print(f"Omega = {omega}")
        print("-" * 78)

        ok3, actual, predicted = check_theorem3(omega)
        all_ok &= ok3
        print(f"  [1] Theorem 3 support characterization: "
              f"{'OK' if ok3 else 'FAILED'}")
        if not ok3:
            print(f"      actual support:    {sorted(actual)}")
            print(f"      predicted support: {sorted(predicted)}")

        ok_cross, s_main, G_main, s_vd, G_vd = check_cross_implementation(omega)
        all_ok &= ok_cross
        print(f"  [2] Cross-check vs verify_dissymmetry_omega.py (sympy): "
              f"{'OK' if ok_cross else 'FAILED'}")
        print(f"      s(x) [main]:  {s_main}")
        print(f"      s(x) [sympy]: {s_vd}")
        print(f"      G(x) [main]:  {G_main}")
        print(f"      G(x) [sympy]: {G_vd}")

        ok_crit, detail = check_critical_pair(omega)
        all_ok &= ok_crit
        if detail[0] == 'theorem5':
            _, rho, tau, tau_exact, rel_diff = detail
            print(f"  [3] Theorem 5 closed form vs numeric critical point: "
                  f"{'OK' if ok_crit else 'FAILED'}")
            print(f"      rho = {nstr(rho, 10)}   tau (numeric) = {nstr(tau, 10)}"
                  f"   tau (closed form) = {nstr(tau_exact, 10)}"
                  f"   relative diff = {nstr(rel_diff, 4)}")
        else:
            _, rho, tau, lhs, diff = detail
            print(f"  [3] Proposition 6 identity at numeric critical pair: "
                  f"{'OK' if ok_crit else 'FAILED'}")
            print(f"      rho = {nstr(rho, 10)}   tau = {nstr(tau, 10)}"
                  f"   LHS = {nstr(lhs, 10)} (should be 2)"
                  f"   |LHS-2| = {nstr(diff, 4)}")

    print()
    print("=" * 78)
    print(f"ALL CHECKS PASSED: {all_ok}")
    print("=" * 78)
    return all_ok


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
