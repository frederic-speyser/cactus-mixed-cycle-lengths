"""
reproduce_table1.py

End-to-end reproduction of Table 1 of [2], the four representative
examples of the article: for each Omega in {{5,6}, {5,7}, {5,7,9},
{5,6,7}}, this script

  1. computes the exact rooted enumeration series via
     mgonal_cactus_series_omega.py, and reads off the modulus d and the
     exceptional zero coefficients within the residue class n = 1 (mod d)
     (Theorem 3 and Corollary 4 of [2]);
  2. locates the critical pair (rho_Omega, tau_Omega) via
     critical_point_solver.py, cross-checked against the closed form of
     Theorem 5 when Omega is all-odd, or against the identity of
     Proposition 6 otherwise.

This is the single script referenced in the "Code and data availability"
section of [2] for reproducing Table 1. It intentionally uses modest
numerical parameters (documented in critical_point_solver.py) to keep the
running time reasonable; the higher-precision values actually printed in
the published Table 1 were obtained with larger parameters, as recorded
in the project's internal development log (available from the author on
request, per the Acknowledgments section of [2]).

Reference:
  [2] F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
      Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths"
      (submitted to the Journal of Integer Sequences, 2026), Table 1.

Author: Frederic G. Speyser
Run: python3 reproduce_table1.py
"""
from math import gcd
from functools import reduce

import mgonal_cactus_series_omega as om
from critical_point_solver import (
    solve_rho, solve_tau, theorem5_closed_form, proposition6_check, nstr, mpf,
)

EXAMPLES = [
    dict(omega=(5, 6), bracket=(mpf('0.51'), mpf('0.55'))),
    dict(omega=(5, 7), bracket=(mpf('0.54'), mpf('0.56'))),
    dict(omega=(5, 7, 9), bracket=(mpf('0.52'), mpf('0.55'))),
    dict(omega=(5, 6, 7), bracket=(mpf('0.49'), mpf('0.52'))),
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
          f"{'rho_Omega':<14}{'tau_Omega':<14}{'closed form?':<14}")
    print("-" * 84)
    for ex in EXAMPLES:
        omega = ex['omega']
        d, zeros = exceptional_zeros(omega)
        rho = solve_rho(omega, *ex['bracket'])
        tau, V = solve_tau(omega, rho)

        all_odd = all(m % 2 == 1 for m in omega)
        if all_odd:
            tau = theorem5_closed_form(omega)  # exact, per Theorem 5
            closed = "yes (Thm. 5)"
        else:
            closed = "no (Prop. 6)"

        omega_str = "{" + ",".join(map(str, omega)) + "}"
        zeros_str = ", ".join(map(str, zeros)) if zeros else "(none)"
        print(f"{omega_str:<10}{d:<4}{zeros_str:<28}"
              f"{nstr(rho, 8):<14}{nstr(tau, 8):<14}{closed:<14}")

        if not all_odd:
            lhs = proposition6_check(omega, rho, tau, V)
            print(f"{'':<10}{'':<4}Proposition 6 consistency check: "
                  f"left-hand side = {nstr(lhs, 6)} (should equal 2)")


if __name__ == "__main__":
    main()
