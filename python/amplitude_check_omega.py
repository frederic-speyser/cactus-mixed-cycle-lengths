"""
Numerical verification of the amplitude-constant identity underlying the
asymptotic coefficient estimate s_n ~ C_Omega * rho_Omega^{-n} * n^{-3/2}
for the rooted enumeration series s_Omega(x).

The amplitude constant C_Omega is computed independently via two routes
and the two values are compared:

  C_dir = d * beta / (2 * sqrt(pi)),
      beta = sqrt(2 * rho * F_x(rho, tau) / F_yy(rho, tau))
      -- original coordinates, F(x,y) = x * exp(K_Omega(x,y) + H_Omega(x))

  C_red = gamma * rho * d**1.5 / (2 * sqrt(pi)),
      gamma = sqrt(2 * z_c * Psi_z / Psi_HH)
      -- reduced coordinates, z = x^d, h(z) = Psi(z, h(z))

F_x, F_yy and Psi_z are evaluated by high-precision Taylor-series
differentiation (mpmath.diff), which is well conditioned at the working
precision used here; a naive fixed-step finite-difference formula loses
accuracy on the second derivative F_yy because the step is squared in
the denominator. Psi_HH is evaluated from its closed form, which holds
for both parities of Omega:

    Psi_HH(z_c, h_c) = rho/tau + (rho/2) * sum_m (m-1)(m-2) tau**(m-2)

This script reuses critical_point_solver.py (high-precision solver for
the critical pair (rho_Omega, tau_Omega)) without duplicating its logic.

Authors: Frederic G. Speyser & Kseniya Vyatkina
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina. All rights reserved.

Run: python3 amplitude_check_omega.py
"""
import math
from math import gcd
from functools import reduce
import mpmath
from mpmath import mp, mpf, nstr, sqrt, pi

from critical_point_solver import (
    compute_series, series_eval, solve_critical_pair, _K_xy, _H,
)

mp.dps = 50


def _d(omega):
    return reduce(gcd, [m - 1 for m in omega])


def _F(omega, s, N, x, y):
    """F(x,y) = x * exp(K_Omega(x,y) + H_Omega(x))."""
    return x * mpmath.e ** (_K_xy(omega, s, N, x, y) + _H(omega, s, N, x))


def _h_coeffs(s, d, N):
    """h(z) = s(x)/x with x^d = z: h_k = s_{1+d k}, k = 0..floor((N-1)/d).
    Note: h_0 = s_1 = 1 is a nonzero constant term -- unlike s, whose
    index 0 is unused/zero. series_eval() (designed for s) skips index 0
    and must not be reused for h; h_eval() below is used instead."""
    Kmax = (N - 1) // d
    h = [mpf(0)] * (Kmax + 1)
    for k in range(Kmax + 1):
        n = 1 + d * k
        if n <= N:
            h[k] = s[n]
    return h


def h_eval(h, Kmax, z):
    """Truncated sum_{k=0}^Kmax h_k z^k, including the constant term h_0."""
    tot = mpf(0)
    zp = mpf(1)
    for k in range(Kmax + 1):
        tot += h[k] * zp
        zp *= z
    return tot


def _Psi(omega, s, h, d, N, Nh, z, H):
    """Psi(z,H) = exp(G_1(z,H) + R(z)); R(z) is H_Omega(x) evaluated at
    x = z**(1/d) (same quantity in reduced coordinates: reuses _H)."""
    zc2 = z * z
    hz2 = h_eval(h, Nh, zc2)  # h(z^2), subcritical, well conditioned
    tot = mpf(0)
    for m in omega:
        qm = (m - 1) // d
        if m % 2 == 1:
            tot += mpf('0.5') * z ** qm * (H ** (m - 1) + hz2 ** ((m - 1) // 2))
        else:
            tot += mpf('0.5') * z ** qm * (H ** (m - 1) + H * hz2 ** ((m - 2) // 2))
    x_of_z = z ** (mpf(1) / d)
    Rz = _H(omega, s, N, x_of_z)
    return mpmath.e ** (tot + Rz)


def check(omega):
    d = _d(omega)
    guess = {
        (5, 6): ('0.531', '0.712'), (5, 7): ('0.550', '0.727'),
        (5, 7, 9): ('0.535', '0.692'), (5, 6, 7): ('0.506', '0.665'),
    }[omega]
    rho, tau, s, N, r1, r2 = solve_critical_pair(omega, *guess)

    # --- C_dir: original coordinates ---
    # mpmath.diff (Taylor/Richardson differentiation) is used instead of a
    # manual finite-difference formula: a fixed step of 1e-25 causes
    # catastrophic cancellation on the SECOND derivative (step squared in
    # the denominator, at the 50-significant-digit working precision).
    Fx = mpmath.diff(lambda x: _F(omega, s, N, x, tau), rho)
    Fyy = mpmath.diff(lambda y: _F(omega, s, N, rho, y), tau, 2)
    beta = sqrt(2 * rho * Fx / Fyy)
    C_dir = d * beta / (2 * sqrt(pi))

    # --- C_red: reduced coordinates ---
    z_c = rho ** d
    h_c = tau / rho
    Nh = (N - 1) // d
    h = _h_coeffs(s, d, N)
    Psi_z = mpmath.diff(lambda z: _Psi(omega, s, h, d, N, Nh, z, h_c), z_c)
    Psi_HH = rho / tau + (rho / 2) * sum((m - 1) * (m - 2) * tau ** (m - 2) for m in omega)
    gamma = sqrt(2 * z_c * Psi_z / Psi_HH)
    C_red = gamma * rho * mpf(d) ** mpf('1.5') / (2 * sqrt(pi))

    return dict(omega=omega, d=d, rho=rho, tau=tau, r1=r1, r2=r2,
                Fx=Fx, Fyy=Fyy, beta=beta, C_dir=C_dir,
                z_c=z_c, h_c=h_c, Psi_z=Psi_z, Psi_HH=Psi_HH, gamma=gamma,
                C_red=C_red, diff=abs(C_red - C_dir))


if __name__ == "__main__":
    for omega in [(5, 6), (5, 7), (5, 7, 9), (5, 6, 7)]:
        r = check(omega)
        print(f"Omega={r['omega']}  d={r['d']}")
        print(f"  rho={nstr(r['rho'],15)}  tau={nstr(r['tau'],15)}  "
              f"(residuals F1={nstr(r['r1'],3)}, F2={nstr(r['r2'],3)})")
        print(f"  F_x(rho,tau)={nstr(r['Fx'],12)}   F_yy(rho,tau)={nstr(r['Fyy'],12)}   "
              f"beta={nstr(r['beta'],12)}   C_dir = d*beta/(2 sqrt(pi)) = {nstr(r['C_dir'],15)}")
        print(f"  z_c={nstr(r['z_c'],12)}  h_c={nstr(r['h_c'],12)}  "
              f"Psi_z={nstr(r['Psi_z'],12)}  Psi_HH={nstr(r['Psi_HH'],12)}  "
              f"gamma={nstr(r['gamma'],12)}")
        print(f"  C_red = gamma*rho*d^1.5/(2 sqrt(pi)) = {nstr(r['C_red'],15)}")
        print(f"  |C_red - C_dir| = {nstr(r['diff'],6)}")
        print()
