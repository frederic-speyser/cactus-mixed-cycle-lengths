"""
mgonal_cactus_series_omega.py

Computes the rooted and unrooted enumeration series for strict cactus
graphs admitting a finite set Omega of cycle lengths (free / non-plane
case), by direct generalization of mgonal_cactus_series.py from [1].

For a singleton Omega = {m}, this reduces exactly to the m-gonal case of
[1]: the kernel K_C becomes a sum of one USEQ term per size in Omega,
instead of a single term (Section 5.1 of [1], generalized). The rest of
the machinery (the outer MSET exponential, the fixed-point equation for
s(x), the dissymmetry theorem for the unrooted series) is structurally
unchanged -- only K_C, and the re-rooted term T_Cm in the dissymmetry
decomposition, become sums over m in Omega instead of single terms.

The correctness of this generalization, and of the exact support
characterization it makes possible, is proved in [2] (Lemma 2, Theorem 3,
Corollary 4). This script's output is used to produce Table 1 of [2]; see
reproduce_table1.py in this directory for the end-to-end reproduction.

Reference:
  [1] F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
      Non-Plane m-Gonal Cactus Graphs via Split-Decomposition", Sections
      5.1-5.3. This script generalizes mgonal_cactus_series.py from the
      repository accompanying that paper.
  [2] F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
      Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths"
      (submitted to the Journal of Integer Sequences, 2026).

Author: Frederic G. Speyser
Run: python3 mgonal_cactus_series_omega.py --omega 5,6
"""
import argparse
from fractions import Fraction as F

N = 100  # truncation order; raised automatically if needed to reach 25 terms


# ---------- formal power series utilities (lists of Fraction, index = degree) ----------

def zero():
    return [F(0)] * (N + 1)


def mul(a, b):
    c = zero()
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        maxj = N - i
        if maxj < 0:
            continue
        for j, bj in enumerate(b[:maxj + 1]):
            if bj == 0:
                continue
            c[i + j] += ai * bj
    return c


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def scale(a, k):
    return [x * k for x in a]


def stretch(a, r):
    """coefficients of a(x^r)"""
    c = zero()
    for n, an in enumerate(a):
        if n * r <= N:
            c[n * r] = an
    return c


def shift_by_x(a):
    """multiply by x"""
    c = zero()
    for n in range(N):
        c[n + 1] = a[n]
    return c


def power_int(a, k):
    r = [F(0)] * (N + 1)
    r[0] = F(1)
    base = a
    while k > 0:
        if k & 1:
            r = mul(r, base)
        base = mul(base, base)
        k >>= 1
    return r


def exp_series(u):
    """exp(u) for a series u with u[0] == 0, via v' = u' v (Cauchy product recurrence)."""
    assert u[0] == 0
    v = zero()
    v[0] = F(1)
    for n in range(1, N + 1):
        s = F(0)
        for k in range(1, n + 1):
            if u[k] != 0:
                s += k * u[k] * v[n - k]
        v[n] = s / n
    return v


# ---------- combinatorial specification, generalized to a set Omega ----------

def K_C_single(s, m):
    """K_C contribution of a single cycle length m (Eq. 2 of [1], Section 5.1)."""
    s2 = stretch(s, 2)
    if m % 2 == 1:
        term1 = power_int(s, m - 1)
        term2 = power_int(s2, (m - 1) // 2)
        return scale(add(term1, term2), F(1, 2))
    else:
        term1 = power_int(s, m - 1)
        term2 = mul(s, power_int(s2, (m - 2) // 2))
        return scale(add(term1, term2), F(1, 2))


def K_C(s, omega):
    """K_C(x) for a mixed Omega: a sum of one USEQ term per size in Omega."""
    total = zero()
    for m in omega:
        total = add(total, K_C_single(s, m))
    return total


def sum_i_KC_xi_over_i(s, omega):
    """Sigma_{i>=1} K_C(x^i)/i, truncated."""
    total = zero()
    min_deg = min(m - 1 for m in omega)
    i = 1
    while i * min_deg <= N:
        s_xi = stretch(s, i)
        kc_i = K_C(s_xi, omega)
        total = add(total, scale(kc_i, F(1, i)))
        i += 1
    return total


def solve_s(omega, iters=None):
    """Fixed point s = x + x*(exp(Sigma K_C(x^i)/i) - 1)."""
    if iters is None:
        iters = N + 2
    s = zero()
    s[1] = F(1)  # start from s = x
    for _ in range(iters):
        E = exp_series(sum_i_KC_xi_over_i(s, omega))
        Eminus1 = list(E)
        Eminus1[0] -= 1
        s_new = zero()
        s_new[1] += 1
        s_new = add(s_new, shift_by_x(Eminus1))
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


def Z_Dm(s, m):
    """Z_{D_m}(s(x), s(x^2), ..., s(x^m)) -- dihedral cycle index (Eq. 8 of [1])."""
    p = {i: stretch(s, i) for i in range(1, m + 1)}
    total = zero()
    for d in divisors(m):
        term = power_int(p[d], m // d)
        total = add(total, scale(term, F(phi_totient(d), 2 * m)))
    if m % 2 == 1:
        extra = mul(p[1], power_int(p[2], (m - 1) // 2))
        total = add(total, scale(extra, F(1, 2)))
    else:
        extra1 = mul(power_int(p[1], 2), power_int(p[2], (m - 2) // 2))
        extra2 = power_int(p[2], m // 2)
        total = add(total, scale(add(extra1, extra2), F(1, 4)))
    return total


def solve_G(omega):
    """Unrooted series via the dissymmetry theorem, generalized to mixed Omega:
    T_Cm becomes a SUM of dihedral cycle indices, one per size in Omega."""
    s = solve_s(omega)
    KC = K_C(s, omega)
    E = exp_series(sum_i_KC_xi_over_i(s, omega))
    Eminus1 = list(E)
    Eminus1[0] -= 1
    S_X = shift_by_x(Eminus1)                    # S_X(x) = x*(E(x)-1)
    S_C = sub(sub(E, [F(1)] + [F(0)] * N), KC)     # S_C(x) = E(x) - 1 - K_C(x)
    T_S = shift_by_x(S_C)                          # T_S(x) = x * S_C(x)
    T_SCm = mul(KC, S_X)                            # T_{S-Cm}(x) = K_C(x) * S_X(x)
    T_Cm = zero()
    for m in omega:
        T_Cm = add(T_Cm, Z_Dm(s, m))
    G = add(sub(T_Cm, T_SCm), T_S)
    return s, G


def series_terms(coeffs, count):
    """First `count` non-zero terms as (degree, int coefficient) pairs.
    Raises if any non-zero coefficient is not a non-negative integer."""
    terms = []
    for n, c in enumerate(coeffs):
        if c == 0:
            continue
        assert c.denominator == 1, f"non-integer coefficient at x^{n}: {c}"
        assert c >= 0, f"negative coefficient at x^{n}: {c}"
        terms.append((n, c.numerator))
        if len(terms) == count:
            break
    return terms


def format_terms(terms):
    parts = []
    for n, c in terms:
        if n == 0:
            parts.append(f"{c}")
        elif n == 1:
            parts.append(f"{c}x" if c != 1 else "x")
        else:
            parts.append(f"{c}x^{n}" if c != 1 else f"x^{n}")
    return " + ".join(parts)


def main():
    global N
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=str, default="5,6",
                         help="Comma-separated cycle lengths, e.g. 5,6")
    parser.add_argument("--terms", type=int, default=25,
                         help="Number of non-zero terms to report")
    args = parser.parse_args()
    omega = tuple(sorted(int(x) for x in args.omega.split(",")))
    assert all(m >= 5 for m in omega), "This script covers m >= 5 only (Lemma 1 of [1])."

    want = args.terms
    # Enlarge N until we have enough non-zero terms in both series.
    while True:
        s, G = solve_G(omega)
        rooted_terms = series_terms(s, want)
        unrooted_terms = series_terms(G, want)
        if len(rooted_terms) >= want and len(unrooted_terms) >= want:
            break
        N *= 2

    print(f"Omega = {{{', '.join(map(str, omega))}}}")
    print(f"Truncation order used: N = {N}\n")
    print(f"rooted (offset 0)  : {format_terms(rooted_terms)}")
    print()
    print(f"unrooted (offset 1): {format_terms(unrooted_terms)}")

    print("\nSanity check: all reported coefficients are non-negative integers. OK.")


if __name__ == "__main__":
    main()
