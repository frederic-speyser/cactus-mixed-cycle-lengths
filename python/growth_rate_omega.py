"""
growth_rate_omega.py

Estimates rho_Omega directly from the exact rational coefficients
computed by mgonal_cactus_series_omega.py, via a coefficient-ratio test
-- a method entirely independent of critical_point_omega.py, which
instead solves the criticality condition Phi_y(rho,tau)=1 numerically by
bisection on a double-precision series evaluation.

The two methods share no code path: this one works from exact integer
coefficients (fractions.Fraction, via mgonal_cactus_series_omega.py) and
a plain ratio of consecutive non-zero terms; critical_point_omega.py
works from a double-precision truncated series and an analytic
criticality condition. Agreement between the two is a genuine
cross-check, not a restatement of the same computation.

Caveat, inherited from the pure-m case and worth restating here: for a
periodic support (single m), the raw ratio test needs many terms to
settle, because of the n^(-3/2) sub-exponential correction. For mixed
Omega the support is aperiodic almost immediately (see the aperiodicity
argument sketched in [1]'s conclusion, and confirmed computationally in
this repository's first entries), so consecutive non-zero terms are
available from very early on -- but the same n^(-3/2) correction still
means the raw ratio converges to 1/rho_Omega only slowly, from above.

Reference: Fr. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
Non-Plane m-Gonal Cactus Graphs via Split-Decomposition" [1], Theorems
2-4 (growth rate); mgonal_cactus_growth_rate.py from that paper's
repository (the pure-m version of this script).

Authors: Frederic G. Speyser & Kseniya Vyatkina
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina. All rights reserved.
Run: python3 growth_rate_omega.py --omega 5,6
"""
import argparse
import time
import mgonal_cactus_series_omega as solver
from critical_point_solver import solve_critical_pair, nstr

# Rough initial guesses for Newton's method, one or two correct digits
# suffice; used only as a starting point for --omega values not in this
# table (the four examples of Table 1 of [2]).
_GUESS = {
    (5, 6): ('0.531', '0.712'),
    (5, 7): ('0.550', '0.727'),
    (5, 7, 9): ('0.535', '0.692'),
    (5, 6, 7): ('0.506', '0.665'),
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=str, default="5,6")
    parser.add_argument("--N", type=int, default=150,
                         help="Truncation order for the ratio test")
    args = parser.parse_args()
    omega = tuple(sorted(int(x) for x in args.omega.split(",")))

    t0 = time.time()
    solver.N = args.N
    s = solver.solve_s(omega)
    nz = [(n, c) for n, c in enumerate(s) if c != 0]
    print(f"Omega = {{{', '.join(map(str, omega))}}}")
    print(f"Non-zero rooted terms found up to N={solver.N}: {len(nz)}")
    print(f"[{time.time()-t0:.1f}s]\n")

    print("Ratio test 1/rho_Omega ~ (s_{n2}/s_{n1})^(1/(n2-n1)), using the")
    print("LAST FEW consecutive pairs of non-zero terms (best available")
    print("approximation at this truncation order):\n")
    for (n1, c1), (n2, c2) in zip(nz[-6:-1], nz[-5:]):
        gap = n2 - n1
        rho_est = (float(c1) / float(c2)) ** (1.0 / gap)
        print(f"  n={n1:4d} -> n={n2:4d}  (gap={gap})   "
              f"rho_Omega (ratio est.) = {rho_est:.6f}   "
              f"1/rho = {1/rho_est:.6f}")

    print(f"\nFor comparison, critical_point_solver.py gives (independent")
    print(f"method: high-precision series-plus-Newton solution of the")
    print(f"characteristic system):")
    guess = _GUESS.get(omega, ('0.5', '0.6'))
    rho, tau, *_ = solve_critical_pair(omega, *guess)
    print(f"  rho_Omega = {nstr(rho, 10)}   1/rho_Omega = {nstr(1/rho, 10)}")

    print(f"\nPure cases from [1], Table 3, for reference:")
    print(f"  rho_5 = 0.604765   1/rho_5 = {1/0.604765:.6f}")
    print(f"  rho_6 = 0.633235   1/rho_6 = {1/0.633235:.6f}")


if __name__ == "__main__":
    main()
