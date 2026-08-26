"""
verify_dissymmetry_omega.py

Independent verification of the unrooted series G(x) for mixed Omega, via
the dissymmetry theorem (Section 5.3 of [1]), reimplemented from scratch
using sympy.Rational arithmetic and hand-written convolutions -- a
different code path from mgonal_cactus_series_omega.py, which uses
fractions.Fraction. Only the idea (fixed-point equation, dissymmetry
decomposition) is shared; every function below is a fresh implementation,
not an import.

An earlier version of this script attempted a fully symbolic SymPy
approach (building and repeatedly sp.expand()-ing polynomial expressions
in x) and was abandoned as impractically slow -- documented here because
it is a genuine methodological choice, not an oversight: manual
convolution with sympy.Rational scalars is both a legitimate independent
implementation (distinct arithmetic type from Fraction, distinct code)
and fast enough to run at the truncation orders needed.

G(x) = T_Cm(x) + T_S(x) - T_{S-Cm}(x)      (dissymmetry theorem, Eq. 4 of [1])
  T_Cm(x)    = sum over m in Omega of Z_Dm(s(x), s(x^2), ..., s(x^m))
  T_S(x)     = x * S_C(x)
  T_{S-Cm}(x) = K_C(x) * S_X(x)

Reference: Fr. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
Non-Plane m-Gonal Cactus Graphs via Split-Decomposition" [1], Section 5.3.

Author: Frederic G. Speyser
Run: python3 verify_dissymmetry_omega.py   (requires: pip install sympy)
"""
import argparse
from sympy import Rational, Integer


def zero(N):
    return [Integer(0)] * (N + 1)


def mul(a, b, N):
    c = zero(N)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        maxj = N - i
        if maxj < 0:
            continue
        for j in range(maxj + 1):
            bj = b[j]
            if bj == 0:
                continue
            c[i + j] += ai * bj
    return c


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def scale(a, k):
    return [v * k for v in a]


def stretch(a, r, N):
    c = zero(N)
    for n, an in enumerate(a):
        if n * r <= N:
            c[n * r] = an
    return c


def power_int(a, k, N):
    r = zero(N)
    r[0] = Integer(1)
    base = a
    while k > 0:
        if k & 1:
            r = mul(r, base, N)
        base = mul(base, base, N)
        k >>= 1
    return r


def exp_series(u, N):
    """exp(u), u[0]=0, via the standard n*v_n = sum k*u_k*v_{n-k} recurrence."""
    v = zero(N)
    v[0] = Integer(1)
    for n in range(1, N + 1):
        acc = Integer(0)
        for k in range(1, n + 1):
            if u[k] != 0:
                acc += k * u[k] * v[n - k]
        v[n] = acc / n
    return v


def K_C_single(s, m, N):
    s2 = stretch(s, 2, N)
    if m % 2 == 1:
        term1 = power_int(s, m - 1, N)
        term2 = power_int(s2, (m - 1) // 2, N)
        return scale(add(term1, term2), Rational(1, 2))
    else:
        term1 = power_int(s, m - 1, N)
        term2 = mul(s, power_int(s2, (m - 2) // 2, N), N)
        return scale(add(term1, term2), Rational(1, 2))


def K_C(s, omega, N):
    total = zero(N)
    for m in omega:
        total = add(total, K_C_single(s, m, N))
    return total


def sum_i_KC_xi_over_i(s, omega, N):
    total = zero(N)
    min_deg = min(m - 1 for m in omega)
    i = 1
    while i * min_deg <= N:
        s_xi = stretch(s, i, N)
        kc_i = K_C(s_xi, omega, N)
        total = add(total, scale(kc_i, Rational(1, i)))
        i += 1
    return total


def solve_s(omega, N, iters=None):
    if iters is None:
        iters = N + 2
    s = zero(N)
    s[1] = Integer(1)
    for _ in range(iters):
        E = exp_series(sum_i_KC_xi_over_i(s, omega, N), N)
        s_new = zero(N)
        s_new[1] = Integer(1)
        for n in range(N):
            s_new[n + 1] += E[n] - (Integer(1) if n == 0 else Integer(0))
        s = s_new
    return s


def phi_totient(n):
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


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def Z_Dm(s, m, N):
    p = {i: stretch(s, i, N) for i in range(1, m + 1)}
    total = zero(N)
    for d in divisors(m):
        term = power_int(p[d], m // d, N)
        total = add(total, scale(term, Rational(phi_totient(d), 2 * m)))
    if m % 2 == 1:
        extra = mul(p[1], power_int(p[2], (m - 1) // 2, N), N)
        total = add(total, scale(extra, Rational(1, 2)))
    else:
        extra1 = mul(power_int(p[1], 2, N), power_int(p[2], (m - 2) // 2, N), N)
        extra2 = power_int(p[2], m // 2, N)
        total = add(total, scale(add(extra1, extra2), Rational(1, 4)))
    return total


def solve_G(omega, N):
    s = solve_s(omega, N)
    KC = K_C(s, omega, N)
    E = exp_series(sum_i_KC_xi_over_i(s, omega, N), N)
    Eminus1 = list(E)
    Eminus1[0] -= 1
    S_X = zero(N)
    for n in range(N):
        S_X[n + 1] = Eminus1[n]
    S_C = sub(Eminus1, KC)
    T_S = zero(N)
    for n in range(N):
        T_S[n + 1] = S_C[n]
    T_SCm = mul(KC, S_X, N)
    T_Cm = zero(N)
    for m in omega:
        T_Cm = add(T_Cm, Z_Dm(s, m, N))
    G = add(sub(T_Cm, T_SCm), T_S)
    return s, G


def nonzero_terms(coeffs, count):
    out = []
    for n, c in enumerate(coeffs):
        if c == 0:
            continue
        assert c.is_Integer and c >= 0, f"bad coefficient at x^{n}: {c}"
        out.append((n, int(c)))
        if len(out) == count:
            break
    return out


def fmt(terms):
    return " + ".join(
        f"{c}x^{n}" if (c != 1 or n == 0) else f"x^{n}" for n, c in terms)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=str, default="5,6")
    parser.add_argument("--N", type=int, default=32)
    parser.add_argument("--terms", type=int, default=8)
    args = parser.parse_args()
    omega = tuple(sorted(int(v) for v in args.omega.split(",")))
    N = args.N

    print(f"Omega = {{{', '.join(map(str, omega))}}}")
    print(f"Assembling s(x) and G(x) independently (sympy.Rational, "
          f"hand-written convolutions) up to N={N}...")
    s, G = solve_G(omega, N)

    got = nonzero_terms(G, args.terms)
    print(f"\nUnrooted series G(x), first {len(got)} non-zero terms "
          f"(independent SymPy implementation):")
    print("  " + fmt(got))

    print(f"\nCross-check against mgonal_cactus_series_omega.py's unrooted")
    print(f"output for Omega={{5,6}} (Fraction-based convolutions, a")
    print(f"different code path):")
    published = [(5, 1), (6, 1), (9, 1), (10, 1), (11, 1), (13, 3), (14, 6), (15, 6)]
    print("  " + fmt(published[:len(got)]))
    match = got == published[:len(got)]
    print(f"\nMatch: {match}")


if __name__ == "__main__":
    main()
