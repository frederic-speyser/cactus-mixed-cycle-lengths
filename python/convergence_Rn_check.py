"""
Coefficient-level convergence check for the amplitude constant:
verifies directly from the rooted enumeration coefficients s_n that

    R_n := s_n * rho_Omega^n * n^{3/2}  -->  C_Omega        (n -> infinity)

for the four Omega families listed in Table 1, at the admissible
indices n = 1 + d*k (d = gcd{m-1 : m in Omega}).

This check is deliberately independent of the critical-point solver at
the coefficient level: s_n is produced by the combinatorial recurrence
(compute_series, block decomposition) alone, never by evaluating a
series near its own radius of convergence. rho_Omega and C_Omega used
for comparison are the values already established by
critical_point_solver.py and amplitude_check_omega.py; solve_critical_pair
is called here to (a) obtain the coefficients s_n needed for R_n and (b)
verify, by explicit assertion, that it reproduces those already-computed
values.

This is a numerical convergence check, not a proof of the rate of
convergence: no claim is made here about the order of the sub-dominant
correction term.

Run: python3 convergence_Rn_check.py

Authors: Frederic G. Speyser & Kseniya Vyatkina
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina. All rights reserved.
"""
import mpmath
from mpmath import mp, mpf, nstr

from critical_point_solver import compute_series, solve_critical_pair
from amplitude_check_omega import check as amplitude_check, _d

mp.dps = 30  # working precision for this convergence check (lower than
             # the 50-digit precision used for the critical-pair solve)

TARGETS = [50, 100, 150, 200, 300]

# Starting points for Newton's method only (approximate; Newton
# reconverges to full precision from any reasonable starting point).
INITIAL_GUESS = {
    (5, 6): ('0.531', '0.712'), (5, 7): ('0.550', '0.727'),
    (5, 7, 9): ('0.535', '0.692'), (5, 6, 7): ('0.506', '0.665'),
}

# Reference values already established by critical_point_solver.py
# (residuals ~1e-51) and amplitude_check_omega.py (|C_red-C_dir| ~1e-52
# or 0), reproduced here for cross-check and reporting.
REFERENCE_VALUES = {
    (5, 6):     dict(rho='0.530923625274025', C='0.156351357866176'),
    (5, 7):     dict(rho='0.549884575826361', C='0.302088001203884'),
    (5, 7, 9):  dict(rho='0.535218127227685', C='0.269362679524554'),
    (5, 6, 7):  dict(rho='0.506181465494348', C='0.138735005523206'),
}
CONSISTENCY_TOL = mpf('1e-13')  # limited by the 15-significant-digit
                                 # truncation of REFERENCE_VALUES above


def run(omega, Nmax=305):
    d = _d(omega)
    rho, tau, s, N, r1, r2 = solve_critical_pair(omega, *INITIAL_GUESS[omega], N=Nmax)
    amp = amplitude_check(omega)
    C = amp['C_dir']  # equals C_red to the working precision

    ref = REFERENCE_VALUES[omega]
    assert abs(rho - mpf(ref['rho'])) < CONSISTENCY_TOL, \
        f"recomputed rho does not match the reference value for {omega}"
    assert abs(C - mpf(ref['C'])) < CONSISTENCY_TOL, \
        f"recomputed C_Omega does not match the reference value for {omega}"

    rows = []
    for target in TARGETS:
        k = round((target - 1) / d)
        n = 1 + d * k
        if n > N:
            continue
        sn = s[n]
        # The support theorem (aperiodicity of the class 1+dk) is an
        # asymptotic statement: s_n != 0 is guaranteed for n large
        # enough, not necessarily at the very first term of the class.
        # No occurrence of s_n == 0 has been observed at the target n
        # values used here (n >= 49); skip gracefully if it occurred.
        if sn == 0:
            continue
        Rn = sn * rho ** n * mpf(n) ** mpf('1.5')
        rows.append((n, Rn, abs(Rn - C) / C))
    return dict(omega=omega, d=d, rho=rho, C=C, rows=rows)


if __name__ == "__main__":
    for omega in [(5, 6), (5, 7), (5, 7, 9), (5, 6, 7)]:
        r = run(omega)
        print(f"Omega={r['omega']}  d={r['d']}  rho={nstr(r['rho'],12)}  "
              f"C_Omega={nstr(r['C'],12)}  (consistent with reference values)")
        for n, Rn, rel in r['rows']:
            print(f"  n={n:4d}   R_n = s_n rho^n n^1.5 = {nstr(Rn,12)}   "
                  f"|R_n - C|/C = {nstr(rel,6)}")
        print()
