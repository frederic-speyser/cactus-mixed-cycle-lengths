"""
mgonal_cactus_series_omega_blocks.py

Computes the rooted and unrooted enumeration series for strict cactus
graphs admitting a finite set Omega of cycle lengths, indexed by NUMBER
OF BLOCKS k -- the convention used by the existing OEIS entries for this
family: the general arrays A332648 (rooted) and A332649 (unrooted), by
Andrew Howroyd, and the single-length columns for the specific k values
relevant to [1]/[2]: k=5 rooted (A398033) and unrooted (A397250), k=6
rooted and unrooted (A398034, A398035), k=7 rooted (A397210) and
unrooted (A398575), and k=8 rooted (A397546); k=8 unrooted has not yet
been submitted. Checked directly against OEIS, one A-number at a time,
on 2026-08-26.

Direct univariate formulation. The block-indexed series s(v) = Sum_k a(k)
v^k satisfies

    s(v) = exp( Sum_{i>=1} K(v^i) / i ),   K(v) = Sum_{m in Omega} K_m(v),

with, for m odd,
    K_m(v) = (v/2) * ( s(v)^(m-1) + s(v^2)^((m-1)/2) ),
and for m even,
    K_m(v) = (v/2) * ( s(v)^(m-1) + s(v) * s(v^2)^((m-2)/2) ).

This is the same construction used throughout this repository (Equation
2 of the companion papers), evaluated directly in the block-counting
variable v rather than tracking vertex count and block count as two
separate variables and marginalizing afterwards. The two routes are
mathematically equivalent -- summing the bivariate series
sum_n [x^n v^k] s(x,v) over all n is the same power series identity as
evaluating the whole construction at x=1, term by term in v, since for
each fixed k only finitely many values of n contribute -- but the
univariate route avoids tracking the vertex dimension at all, which is
the dimension that made the computation expensive. In practice this
turns a computation that would need on the order of 10^2-10^3 hours to
reach 100 terms into one that takes about a minute (up to ~70s for the
heaviest Omega tested, {5,7,9}, measured directly rather than assumed).

Both routes were independently implemented and cross-checked against
each other -- and against the values already published and validated in
the companion exploratory repository cactus-split-decomp-omega for
Omega={5,6} -- before this version replaced the bivariate one as the
one used to generate this repository's data. See CHANGELOG.md for the
detailed record.

The unrooted series follows the same route: writing K_C(v) for the raw
(unmarked) sum above, the dissymmetry theorem gives

    G(v) = T_Cm(v) + T_S(v) - T_{S-Cm}(v)

with T_S(v) = s(v) - 1 - K_C(v) (shifted by one block for the star
itself), T_{S-Cm}(v) = K_C(v) * (s(v) - 1), and T_Cm(v) the dihedral
cycle index of Equation 8, evaluated at p_i = s(v^i) and shifted by one
block for the cycle itself.

Reference: F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths" (in
preparation for submission to the Journal of Integer Sequences, 2026).

Author: Frederic G. Speyser
Run: python3 mgonal_cactus_series_omega_blocks.py --omega 5,6 --terms 100
"""
import argparse
from fractions import Fraction as F


def _make_solver(K):
    """Builds add/mul/power/stretch/exp primitives truncated at degree K
    in the block-counting variable v, and the two top-level solvers
    (solve_s for the rooted series, solve_G for the unrooted one)."""

    def zero():
        return [F(0)] * (K + 1)

    def add(a, b):
        return [x + y for x, y in zip(a, b)]

    def sub(a, b):
        return [x - y for x, y in zip(a, b)]

    def scale(a, c):
        return [x * c for x in a]

    def mul(a, b):
        c = zero()
        for i, ai in enumerate(a):
            if ai == 0:
                continue
            for j, bj in enumerate(b):
                if i + j > K:
                    break
                if bj == 0:
                    continue
                c[i + j] += ai * bj
        return c

    def power_int(a, p):
        r = zero()
        r[0] = F(1)
        base = a
        while p > 0:
            if p & 1:
                r = mul(r, base)
            base = mul(base, base)
            p >>= 1
        return r

    def stretch(a, k):
        c = zero()
        for n, an in enumerate(a):
            if n * k <= K:
                c[n * k] = an
        return c

    def shift1(a):
        c = zero()
        for n in range(K):
            c[n + 1] = a[n]
        return c

    def exp_series(u):
        v = zero()
        v[0] = F(1)
        for n in range(1, K + 1):
            s = F(0)
            for k in range(1, n + 1):
                if u[k] != 0:
                    s += k * u[k] * v[n - k]
            v[n] = s / n
        return v

    def K_C_term(m, s1, s2):
        """K_m evaluated with s(v)=s1, s(v^2)=s2, WITHOUT the leading v
        (that factor is applied by the caller, at the point of use)."""
        if m % 2 == 1:
            t1 = power_int(s1, m - 1)
            t2 = power_int(s2, (m - 1) // 2)
            return scale(add(t1, t2), F(1, 2))
        else:
            t1 = power_int(s1, m - 1)
            t2 = mul(s1, power_int(s2, (m - 2) // 2))
            return scale(add(t1, t2), F(1, 2))

    def K_raw(s, omega):
        total = zero()
        for m in omega:
            total = add(total, K_C_term(m, s, stretch(s, 2)))
        return total

    def solve_s(omega, iters=None):
        """Rooted block-indexed series s(v); iters defaults to K+5,
        empirically sufficient for full correctness at this truncation
        order (verified against independently-known values)."""
        if iters is None:
            iters = K + 5
        s = zero()
        s[0] = F(1)
        for _ in range(iters):
            total = zero()
            i = 1
            while i <= K:
                si = stretch(s, i)
                Ki = K_raw(si, omega)
                Ki_shift = zero()
                for n in range(K + 1 - i):
                    Ki_shift[n + i] = Ki[n]
                total = add(total, scale(Ki_shift, F(1, i)))
                i += 1
            s = exp_series(total)
        return s

    def _phi_totient(n):
        result = n
        p = 2
        nn = n
        while p * p <= nn:
            if nn % p == 0:
                while nn % p == 0:
                    nn //= p
                result -= result // p
            p += 1
        if nn > 1:
            result -= result // nn
        return result

    def _divisors(n):
        return [d for d in range(1, n + 1) if n % d == 0]

    def Z_Dm(s, m):
        """Dihedral cycle index (Equation 8), shifted by one block for
        the cycle itself."""
        p = {i: stretch(s, i) for i in range(1, m + 1)}
        total = zero()
        for d in _divisors(m):
            term = power_int(p[d], m // d)
            total = add(total, scale(term, F(_phi_totient(d), 2 * m)))
        if m % 2 == 1:
            extra = mul(p[1], power_int(p[2], (m - 1) // 2))
            total = add(total, scale(extra, F(1, 2)))
        else:
            extra1 = mul(power_int(p[1], 2), power_int(p[2], (m - 2) // 2))
            extra2 = power_int(p[2], m // 2)
            total = add(total, scale(add(extra1, extra2), F(1, 4)))
        return shift1(total)

    def solve_G(omega, s):
        """La serie G(v) indexee par blocs et non enracinee s'obtient a partir de
        la serie enracinee s (deja calculee) en utilisant le theoreme de dissymetrie."""
        KC = shift1(K_raw(s, omega))
        S_X = sub(s, [F(1)] + [F(0)] * K)   # s(v) - 1
        S_C = sub(S_X, KC)
        T_S = S_C
        T_SCm = mul(KC, S_X)
        T_Cm = zero()
        for m in omega:
            T_Cm = add(T_Cm, Z_Dm(s, m))
        return add(sub(T_Cm, T_SCm), T_S)

    return solve_s, solve_G


def compute(omega, terms=100):
    """Returns (rooted, unrooted): rooted[n] = a(n) for n=0..terms
    (offset 0), unrooted[i] = a(i+1) for i=0..terms-1 (offset 1)."""
    K = terms + 1
    solve_s, solve_G = _make_solver(K)
    s = solve_s(omega)
    G = solve_G(omega, s)
    rooted = [int(v) for v in s[:terms + 1]]
    unrooted = [int(v) for v in G[1:terms + 1]]
    return rooted, unrooted


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=str, default="5,6",
                         help="Comma-separated cycle lengths, e.g. 5,6")
    parser.add_argument("--terms", type=int, default=100,
                         help="Number of terms to compute (a(0) or a(1) through a(terms))")
    args = parser.parse_args()
    omega = tuple(sorted(int(x) for x in args.omega.split(",")))

    rooted, unrooted = compute(omega, args.terms)
    print(f"Omega = {{{', '.join(map(str, omega))}}}")
    print(f"\nrooted (offset 0), {len(rooted)} terms:")
    print(", ".join(map(str, rooted)))
    print(f"\nunrooted (offset 1), {len(unrooted)} terms:")
    print(", ".join(map(str, unrooted)))
