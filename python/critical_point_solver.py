"""
critical_point_solver.py

High-precision, direct numerical solver for the critical pair (rho_Omega,
tau_Omega) of the rooted enumeration series s_Omega(x) of strict cactus
graphs with a finite set Omega of admissible cycle lengths (free /
non-plane case), by a damped fixed-point iteration on the coupled system
s_Omega(x), s_Omega(x^2), s_Omega(x^3), ..., with Aitken extrapolation to
accelerate convergence near the singularity, followed by a square-root
singularity fit to extrapolate tau_Omega = lim_{x -> rho_Omega^-} s_Omega(x).

This consolidates and supersedes two earlier, narrower prototypes: an
odd-Omega-only version, and an independent extension to mixed parities.
The version here handles both cases uniformly, since K_C^{(m)} is given
here for either parity of m directly (Equation 2 of [1]).

Two consistency checks are used throughout the accompanying project and
are reproduced in reproduce_table1.py:
  - for Omega all-odd, tau_Omega found here should match the closed-form
    root of Theorem 5 of [2] to the precision requested;
  - for Omega containing an even length, the critical pair found here
    should satisfy the identity of Proposition 6 of [2] to the precision
    requested.

This script deliberately does not attempt to reach the same precision
used to produce the final published table (obtained with larger K,
more iterations, and more bisection rounds than the defaults below);
see the docstrings of solve_rho and solve_tau for how to trade runtime
for precision.

Reference:
  [1] F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
      Non-Plane m-Gonal Cactus Graphs via Split-Decomposition."
  [2] F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
      Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths"
      (in preparation for submission to the Journal of Integer 
      Sequences, 2026), Theorem 5, Proposition 6, Section 8.

Author: Frederic G. Speyser
Run: python3 critical_point_solver.py
"""
import mpmath
from mpmath import mp, mpf, nstr

mp.dps = 30  # working precision, in significant decimal digits


def K_C_term(m, v1, v2):
    """K_C^{(m)} evaluated at a point where s(t) = v1 and s(t^2) = v2,
    Equation 2 of [1], valid for m of either parity."""
    if m % 2 == 1:
        return (v1 ** (m - 1) + v2 ** ((m - 1) // 2)) / 2
    else:
        return (v1 ** (m - 1) + v1 * v2 ** ((m - 2) // 2)) / 2


def K_omega_at(idx, V, omega, K):
    """K_Omega(x^idx), using the current estimates V[idx] ~ s(x^idx) and
    V[2*idx] ~ s(x^(2*idx))."""
    v1 = V[idx]
    v2 = V[2 * idx] if 2 * idx <= K else mpf(0)
    return sum(K_C_term(m, v1, v2) for m in omega)


def picard_solve_aitken(x, omega, K=50, iters=400, damping=mpf('0.6'), tol=None):
    """Damped fixed-point iteration for s(x^i), i = 1..K, at fixed x,
    with Aitken extrapolation on the i=1 subsequence. Returns
    (V, iterations_used, aitken_estimate_of_s(x)), or (None, it, None) if
    the iteration diverges (indicating x > rho_Omega)."""
    powers = [None] + [x ** i for i in range(1, K + 1)]
    V = [mpf(0)] + [powers[i] for i in range(1, K + 1)]
    tol = tol or mpf(10) ** (-mp.dps + 8)
    hist = []
    prev_aitken = None
    for it in range(iters):
        newV = [mpf(0)] * (K + 1)
        for i in range(1, K + 1):
            if powers[i] < mpf(10) ** (-mp.dps - 5):
                newV[i] = powers[i]
                continue
            expo = mpf(0)
            j = 1
            while i * j <= K:
                expo += K_omega_at(i * j, V, omega, K) / j
                j += 1
            newV[i] = powers[i] * mp.e ** expo
        V = [mpf(0)] + [damping * newV[i] + (1 - damping) * V[i] for i in range(1, K + 1)]
        if not mp.isfinite(V[1]) or V[1] > 50:
            return None, it, None
        hist.append(V[1])
        if len(hist) >= 3:
            a, b, c = hist[-3], hist[-2], hist[-1]
            denom = c - 2 * b + a
            aitken = c - (c - b) ** 2 / denom if denom != 0 else c
            if len(hist) >= 4 and prev_aitken is not None and abs(aitken - prev_aitken) < tol:
                return V, it, aitken
            prev_aitken = aitken
    a, b, c = hist[-3], hist[-2], hist[-1]
    denom = c - 2 * b + a
    aitken = c - (c - b) ** 2 / denom if denom != 0 else c
    return V, iters, aitken


def solve_rho(omega, lo, hi, iters_probe=500, K=40, rounds=22):
    """Locate rho_Omega by bisection: for x < rho_Omega the fixed-point
    iteration converges and stays bounded; for x > rho_Omega it diverges.
    Increase `rounds` for more significant digits (each round halves the
    bracket), and `iters_probe` if convergence is not reached within it
    for x close to rho_Omega."""
    for _ in range(rounds):
        mid = (lo + hi) / 2
        V, it, aitken = picard_solve_aitken(mid, omega, K=K, iters=iters_probe)
        if V is None or aitken is None or aitken > 30 or it >= iters_probe:
            hi = mid
        else:
            lo = mid
    return lo


def solve_tau(omega, rho, epsilons=('0.004', '0.0025', '0.0015', '0.0008'), K=50, iters=4000):
    """Extrapolate tau_Omega = lim_{x -> rho^-} s_Omega(x) by a
    least-squares fit of s(x) = tau - c*sqrt(rho - x) (the square-root
    singularity expected from the smooth implicit-function schema) to
    several points just below rho. Returns (tau, last_V) where last_V is
    the state vector at the closest evaluated point, useful for
    evaluating s(rho^2) as required by Proposition 6's check.
    """
    xs = [rho * (1 - mpf(eps)) for eps in epsilons]
    pts = []
    for x in xs:
        V, it, aitken = picard_solve_aitken(x, omega, K=K, iters=iters)
        pts.append((x, aitken, V))
    S = [mpmath.sqrt(rho - x) for x, _, _ in pts]
    Y = [y for _, y, _ in pts]
    n = len(pts)
    sS = sum(S)
    sY = sum(Y)
    sSS = sum(s * s for s in S)
    sSY = sum(s * y for s, y in zip(S, Y))
    A = mpmath.matrix([[n, -sS], [-sS, sSS]])
    b = mpmath.matrix([sY, -sSY])
    tau, _ = mpmath.lu_solve(A, b)
    return tau, pts[-1][2]


def theorem5_closed_form(omega):
    """Closed-form tau_Omega for Omega all-odd (Theorem 5 of [2]): the
    unique positive root of Sum_{m in Omega} (m-1) y^(m-1) = 2."""
    assert all(m % 2 == 1 for m in omega), "Theorem 5 requires all-odd Omega"
    def poly(y):
        return sum((m - 1) * y ** (m - 1) for m in omega) - 2
    return mpmath.findroot(poly, mpf('0.7'))


def proposition6_check(omega, rho, tau, V_at_rho2):
    """Left-hand side of the identity of Proposition 6 of [2], which
    should equal 2 at the true critical pair (checked only for Omega
    containing an even length; use as a numerical consistency check,
    not as a substitute for the closed form of Theorem 5)."""
    A_rho = sum((V_at_rho2[2] ** ((m - 2) // 2)) / 2 for m in omega if m % 2 == 0)
    lhs = sum((m - 1) * tau ** (m - 1) for m in omega) + 2 * tau * A_rho
    return lhs


if __name__ == "__main__":
    for omega, bracket in [((5, 6), (mpf('0.51'), mpf('0.55'))),
                            ((5, 7), (mpf('0.54'), mpf('0.56'))),
                            ((5, 7, 9), (mpf('0.52'), mpf('0.55'))),
                            ((5, 6, 7), (mpf('0.49'), mpf('0.52')))]:
        rho = solve_rho(omega, *bracket)
        tau, V = solve_tau(omega, rho)
        print(f"Omega={omega}:  rho = {nstr(rho, 10)}   tau = {nstr(tau, 10)}")
        if all(m % 2 == 1 for m in omega):
            tau_exact = theorem5_closed_form(omega)
            print(f"    Theorem 5 closed form: tau = {nstr(tau_exact, 15)}"
                  f"  (relative difference from numeric: {nstr(abs(tau - tau_exact) / tau_exact, 6)})")
        else:
            lhs = proposition6_check(omega, rho, tau, V)
            print(f"    Proposition 6 check: left-hand side = {nstr(lhs, 8)} (should equal 2)")
        print()
