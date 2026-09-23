"""
reproduce_table1.py

Authors: Frederic G. Speyser & Kseniya Vyatkina
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina. All rights reserved.

End-to-end reproduction of Table 1 of [2], the four representative
examples of the article: for each Omega in {{5,6}, {5,7}, {5,7,9},
{5,6,7}}, this script

  1. computes the exact rooted enumeration series via
     mgonal_cactus_series_omega.py, and reads off the modulus d and the
     exceptional zero coefficients within the residue class n = 1 (mod d)
     (Theorem 3 and Corollary 4 of [2]);
  2. locates the critical pair (rho_Omega, tau_Omega) via
     critical_point_solver.py's high-precision series-plus-Newton solver,
     cross-checked against the closed form of Theorem 5 when Omega is
     all-odd, or against the identity of Proposition 6 otherwise.

This is the single script referenced in the "Code and data availability"
section of [2] for reproducing Table 1. Unlike the previous version of
this script, it reaches full working precision (50 significant digits by
default) directly, with no undocumented parameter tuning: the initial
guesses below only need one or two correct digits, since Newton's method
converges quadratically from there. See CHANGELOG.md for the correction
this replaces (the printed Table 1 values for Omega = {5,6} and
{5,6,7} in the originally submitted manuscript did not satisfy the
Proposition 6 identity to the claimed precision; this script's output
does, and matches the values used in the revised manuscript).

Reference:
  [2] F. G. Speyser and K. Vyatkina, "Enumeration and Asymptotic Analysis
      of Strict Non-Plane Cactus Graphs over a Finite Set of Cycle
      Lengths" (submitted to the Journal of Integer Sequences, 2026),
      Table 1.

Run: python3 reproduce_table1.py
"""
from math import gcd
from functools import reduce

import mgonal_cactus_series_omega as om
from critical_point_solver import (
    solve_critical_pair, theorem5_closed_form, proposition6_check, nstr,
)

EXAMPLES = [
    dict(omega=(5, 6), guess=('0.531', '0.712')),
    dict(omega=(5, 7), guess=('0.550', '0.727')),
    dict(omega=(5, 7, 9), guess=('0.535', '0.692')),
    dict(omega=(5, 6, 7), guess=('0.506', '0.665')),
]


def exceptional_zeros(omega, bound=40):
    """The exceptional n = 1 (mod d) with s_n = 0, per Theorem 3 /
    Corollary 4: computed directly from the exact series, independently
    of the Sylvester--Frobenius formula, as a cross-check of Corollary 4
    rather than an application of it."""
    om.N = bound
    s, _ = om.solve_G(omega)
    d = reduce(gcd, [m - 1 for m in omega])
    zeros = [n for n in range(1, bound + 1)
             if (n - 1) % d == 0 and s[n] == 0]
    return d, zeros


def main():
    print(f"{'Omega':<10}{'d':<4}{'exceptional n (s_n=0)':<28}"
          f"{'rho_Omega':<16}{'tau_Omega':<16}{'closed form?':<14}")
    print("-" * 88)
    for ex in EXAMPLES:
        omega = ex['omega']
        d, zeros = exceptional_zeros(omega)
        rho, tau, s, N, r1, r2 = solve_critical_pair(omega, *ex['guess'])

        all_odd = all(m % 2 == 1 for m in omega)
        if all_odd:
            tau_exact = theorem5_closed_form(omega)
            closed = "yes (Thm. 5)"
        else:
            closed = "no (Prop. 6)"

        omega_str = "{" + ",".join(map(str, omega)) + "}"
        zeros_str = ", ".join(map(str, zeros)) if zeros else "(none)"
        print(f"{omega_str:<10}{d:<4}{zeros_str:<28}"
              f"{nstr(rho, 10):<16}{nstr(tau, 10):<16}{closed:<14}")

        if all_odd:
            print(f"{'':<10}{'':<4}Theorem 5 cross-check: numeric tau agrees with "
                  f"closed form to {nstr(abs(tau - tau_exact), 3)}")
        else:
            lhs = proposition6_check(omega, rho, tau, s, N)
            print(f"{'':<10}{'':<4}Proposition 6 consistency check: "
                  f"left-hand side = {nstr(lhs, 12)} (should equal 2)")


if __name__ == "__main__":
    main()
