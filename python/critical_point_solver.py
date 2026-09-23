"""
critical_point_solver.py

Authors: Frederic G. Speyser & Kseniya Vyatkina
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina. All rights reserved.

High-precision, direct numerical solver for the critical pair (rho_Omega,
tau_Omega) of the rooted enumeration series s_Omega(x) of strict cactus
graphs with a finite set Omega of admissible cycle lengths (free,
non-plane, vertex-rooted case).

Method (replaces the earlier damped-fixed-point / Aitken-extrapolation
engine, which was found not to converge to the stated precision for the
mixed-parity examples -- see CHANGELOG.md):

  1. The coefficients s_1, s_2, ..., s_N of s_Omega(x) are computed once,
     exactly (as arbitrary-precision floats, mpmath), via the recursive
     block decomposition of Equation (1) of [2]. This is independent of
     rho_Omega and uses no numerical root-finding.
  2. For any x with |x| < rho_Omega, s_Omega(x) is well approximated by
     the truncated sum sum_{n=1}^N s_n x^n; convergence is geometric away
     from x = rho_Omega, so this is well-conditioned for the "subcritical"
     evaluations s_Omega(rho^2), s_Omega(rho^3), ... that occur inside the
     kernel K_Omega and its tail H_Omega (Section 2 of [2]).
  3. The critical pair (rho, tau) is the solution of the two-equation
     characteristic system of the smooth implicit-function schema
     (Flajolet-Sedgewick, Analytic Combinatorics, Ch. VII):
         F1(rho, tau) := tau - rho * exp(K_Omega(rho, tau) + H_Omega(rho)) = 0
         F2(rho, tau) := sum_{m in Omega} (m-1) tau^(m-1) + 2*tau*A(rho) - 2 = 0
     solved by 2-D Newton's method in high precision (mpmath, default 50
     significant decimal digits), using only the subcritical truncated
     series from step 2 -- never a series evaluated exactly at its own
     radius of convergence, which is the ill-conditioned computation the
     previous engine attempted implicitly.

For Omega all-odd, F2 alone gives the closed-form polynomial of Theorem 5
of [2]; theorem5_closed_form is kept for an independent cross-check.
For Omega containing an even length, proposition6_check evaluates the
left-hand side of the identity of Proposition 6 of [2] at the returned
(rho, tau); it should equal 2 to the working precision.

Reference:
  [1] F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
      Non-Plane m-Gonal Cactus Graphs via Split-Decomposition."
  [2] F. G. Speyser and K. Vyatkina, "Enumeration and Asymptotic Analysis
      of Strict Non-Plane Cactus Graphs over a Finite Set of Cycle
      Lengths" (submitted to the Journal of Integer Sequences, 2026),
      Theorem 5 (thm:odd), Proposition 6 (prop:even), Section 8.

Run: python3 critical_point_solver.py
"""
import mpmath
from mpmath import mp, mpf, nstr

mp.dps = 50  # working precision, in significant decimal digits


# ---------------------------------------------------------------------------
# Step 1-2: the rooted enumeration series, computed once per Omega.
# ---------------------------------------------------------------------------

def compute_series(omega, N=150):
    """Coefficients s[1..N] of s_Omega(x) = sum_n s_n x^n, as mpmath floats.
    s[0] is unused (kept 0) so that s[n] is the coefficient of x^n."""

    def poly_mul_trunc(a, b, deg):
        res = [mpf(0)] * (deg + 1)
        for i, ai in enumerate(a):
            if ai == 0 or i > deg:
                continue
            for j, bj in enumerate(b):
                if bj == 0 or i + j > deg:
                    continue
                res[i + j] += ai * bj
        return res

    def poly_pow_trunc(a, p, deg):
        res = [mpf(0)] * (deg + 1)
        res[0] = mpf(1)
        for _ in range(p):
            res = poly_mul_trunc(res, a, deg)
        return res

    def kappa_upto(d, known_s):
        # K_Omega(x) mod x^{d+1}, using known_s[1..d-1] only (provably
        # sufficient, since every exponent of s(x) in K_C^{(m)} is >= m-1 >= 4).
        sx = known_s[: d + 1]
        sx2 = [mpf(0)] * (d + 1)
        for n in range(1, d + 1):
            if 2 * n <= d:
                sx2[2 * n] = known_s[n]
        tot = [mpf(0)] * (d + 1)
        for m in omega:
            p1 = poly_pow_trunc(sx, m - 1, d)
            if m % 2 == 1:
                e = (m - 1) // 2
                p2 = poly_pow_trunc(sx2, e, d)
                for k in range(d + 1):
                    tot[k] += mpf('0.5') * (p1[k] + p2[k])
            else:
                e = (m - 2) // 2
                p2 = poly_pow_trunc(sx2, e, d)
                conv = poly_mul_trunc(sx, p2, d)
                for k in range(d + 1):
                    tot[k] += mpf('0.5') * (p1[k] + conv[k])
        return tot

    s = [mpf(0)] * (N + 1)
    s[1] = mpf(1)
    for _ in range(N + 2):
        k = kappa_upto(N, s)
        L = [mpf(0)] * (N + 1)
        for i in range(1, N + 1):
            if k[i] == 0:
                continue
            ii = 1
            while i * ii <= N:
                L[i * ii] += k[i] / ii
                ii += 1
        E = [mpf(0)] * (N + 1)
        E[0] = mpf(1)
        for nn in range(1, N + 1):
            acc = mpf(0)
            for j in range(1, nn + 1):
                acc += j * L[j] * E[nn - j]
            E[nn] = acc / nn
        news = [mpf(0)] * (N + 1)
        for nn in range(0, N):
            news[nn + 1] = E[nn]
        diff = max(abs(news[i] - s[i]) for i in range(N + 1))
        s = news
        if diff < mpf(10) ** (-mp.dps + 5):
            break
    return s


def series_eval(s, N, x):
    """Truncated sum_{n=1}^N s_n x^n. Well-conditioned for |x| bounded away
    from the radius of convergence (i.e. for x = rho^i, i >= 2)."""
    tot = mpf(0)
    xp = mpf(1)
    for n in range(1, N + 1):
        xp *= x
        tot += s[n] * xp
    return tot


# ---------------------------------------------------------------------------
# Step 3: the characteristic system.
# ---------------------------------------------------------------------------

def _K_xy(omega, s, N, x, y):
    """Bivariate kernel K_Omega(x, y), Equation (kernel) of [2]: the s(x)
    dependence of the i=1 term of the fixed-point equation is replaced by
    the free variable y; the s(x^2) dependence is evaluated from the
    (subcritical) truncated series."""
    sx2 = series_eval(s, N, x * x)
    tot = mpf(0)
    for m in omega:
        if m % 2 == 1:
            tot += mpf('0.5') * (y ** (m - 1) + sx2 ** ((m - 1) // 2))
        else:
            tot += mpf('0.5') * (y ** (m - 1) + y * sx2 ** ((m - 2) // 2))
    return tot


def _H(omega, s, N, x, imax=60):
    """Tail H_Omega(x) = sum_{i>=2} K_Omega(x^i)/i, using the univariate
    kernel K_Omega(z) = K_Omega(z, s(z)) at z = x^i, i >= 2 -- all
    subcritical when x = rho_Omega < 1, so s(z) is read from the truncated
    series with excellent conditioning."""
    tot = mpf(0)
    for i in range(2, imax + 1):
        z = x ** i
        sz = series_eval(s, N, z)
        term = _K_xy(omega, s, N, z, sz) / i
        tot += term
        if abs(term) < mpf(10) ** (-mp.dps - 5) and i > 5:
            break
    return tot


def _A(omega, s, N, x):
    """A(x) = (1/2) sum_{m in Omega, m even} s(x^2)^((m-2)/2), the extra
    term of Proposition 6 of [2]."""
    sx2 = series_eval(s, N, x * x)
    tot = mpf(0)
    for m in omega:
        if m % 2 == 0:
            tot += mpf('0.5') * sx2 ** ((m - 2) // 2)
    return tot


def solve_critical_pair(omega, guess_rho, guess_tau, N=150, s=None):
    """Solve the characteristic system F1 = F2 = 0 for (rho_Omega,
    tau_Omega) by 2-D Newton's method, starting from (guess_rho,
    guess_tau). Returns (rho, tau, s, N, residual_F1, residual_F2)."""
    if s is None:
        s = compute_series(omega, N)

    def F1(rho, tau):
        return tau - rho * mpmath.e ** (_K_xy(omega, s, N, rho, tau) + _H(omega, s, N, rho))

    def F2(rho, tau):
        tot = sum((m - 1) * tau ** (m - 1) for m in omega)
        tot += 2 * tau * _A(omega, s, N, rho)
        return tot - 2

    rho, tau = mpf(guess_rho), mpf(guess_tau)
    h = mpf(10) ** (-mp.dps // 2)
    for _ in range(60):
        f1, f2 = F1(rho, tau), F2(rho, tau)
        df1dr = (F1(rho + h, tau) - F1(rho - h, tau)) / (2 * h)
        df1dt = (F1(rho, tau + h) - F1(rho, tau - h)) / (2 * h)
        df2dr = (F2(rho + h, tau) - F2(rho - h, tau)) / (2 * h)
        df2dt = (F2(rho, tau + h) - F2(rho, tau - h)) / (2 * h)
        det = df1dr * df2dt - df1dt * df2dr
        drho = -(f1 * df2dt - f2 * df1dt) / det
        dtau = -(df1dr * f2 - df2dr * f1) / det
        rho += drho
        tau += dtau
        if abs(drho) < mpf(10) ** (-mp.dps + 10) and abs(dtau) < mpf(10) ** (-mp.dps + 10):
            break
    return rho, tau, s, N, F1(rho, tau), F2(rho, tau)


def theorem5_closed_form(omega):
    """Closed-form tau_Omega for Omega all-odd (Theorem 5 / thm:odd of
    [2]): the unique positive root of sum_{m in Omega} (m-1) y^(m-1) = 2."""
    assert all(m % 2 == 1 for m in omega), "Theorem 5 requires all-odd Omega"

    def poly(y):
        return sum((m - 1) * y ** (m - 1) for m in omega) - 2

    return mpmath.findroot(poly, mpf('0.7'))


def proposition6_check(omega, rho, tau, s, N):
    """Left-hand side of the identity of Proposition 6 / prop:even of [2],
    which should equal 2 at the true critical pair (checked only for Omega
    containing an even length; a numerical consistency check, not a
    substitute for the closed form of Theorem 5)."""
    lhs = sum((m - 1) * tau ** (m - 1) for m in omega) + 2 * tau * _A(omega, s, N, rho)
    return lhs


if __name__ == "__main__":
    # Initial guesses only need one or two correct digits; Newton's method
    # converges quadratically to the full working precision from there.
    for omega, (g_rho, g_tau) in [
        ((5, 6), ('0.531', '0.712')),
        ((5, 7), ('0.550', '0.727')),
        ((5, 7, 9), ('0.535', '0.692')),
        ((5, 6, 7), ('0.506', '0.665')),
    ]:
        rho, tau, s, N, r1, r2 = solve_critical_pair(omega, g_rho, g_tau)
        print(f"Omega={omega}:  rho = {nstr(rho, 15)}   tau = {nstr(tau, 15)}"
              f"   (residuals: {nstr(r1, 3)}, {nstr(r2, 3)})")
        if all(m % 2 == 1 for m in omega):
            tau_exact = theorem5_closed_form(omega)
            print(f"    Theorem 5 closed form: tau = {nstr(tau_exact, 15)}"
                  f"  (relative difference from numeric: {nstr(abs(tau - tau_exact) / tau_exact, 6)})")
        else:
            lhs = proposition6_check(omega, rho, tau, s, N)
            print(f"    Proposition 6 check: left-hand side = {nstr(lhs, 10)} (should equal 2)")
        print()
