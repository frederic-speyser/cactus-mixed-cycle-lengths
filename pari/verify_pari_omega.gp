/*
verify_pari_omega.gp

A second, fully independent solver for the mixed-Omega case, using PARI's
native truncated power series arithmetic -- a completely different code
path from mgonal_cactus_series_omega.py (Python, fractions.Fraction) and
verify_dissymmetry_omega.py (Python, sympy.Rational): different language,
different runtime, and PARI's own internal representation and truncation
of O(x^N) series rather than any hand-written convolution.

This is the mixed-Omega analogue of verify_pari.gp from [1], generalized
so that K_C becomes a sum of one term per size in Omega (instead of a
single term), and the dissymmetry re-rooted term T_Cm becomes a sum of
dihedral cycle indices, one per size in Omega.

Development notes (kept here because both were genuinely surprising):
  1. "omega" is a PARI/GP built-in (number of distinct prime factors of
     an integer) -- using it as a variable name silently shadows it and
     produces a confusing cascade of unrelated-looking parse errors.
     Renamed to "Om" throughout.
  2. PARI/GP's parser does not support embedded braces: a function body
     already wrapped in {...} cannot contain a further {...}-grouped
     multi-statement for/while loop inside it -- this raises "sorry,
     embedded braces (in parser) is not yet implemented" and, worse,
     corrupts parsing of everything downstream with unrelated-looking
     domain errors. Fixed by rewriting every loop that lived inside a
     function body as a single call to PARI's own sum(), which needs no
     block at all; multi-statement {...} blocks are used only at the
     script's top level, one level deep, never nested inside a function.
  3. gp does not exit after running a script file unless the script
     itself calls quit -- omitting it hangs the process waiting on
     stdin, easily mistaken for an infinite loop in the computation.

Reference: Fr. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
Non-Plane m-Gonal Cactus Graphs via Split-Decomposition" [1], Sections
5.1-5.3.

Author: Frederic G. Speyser
Run: gp -q verify_pari_omega.gp < /dev/null
*/

N = 40;
default(seriesprecision, N + 2);
Om = [5, 6];

Kc(s, m) =
{
  my(s2 = subst(s, x, x^2));
  if(m % 2 == 1,
    1/2 * (s^(m - 1) + s2^((m - 1) / 2))
  ,
    1/2 * (s^(m - 1) + s * s2^((m - 2) / 2))
  );
}

KCsum(s, om) = sum(k = 1, #om, Kc(s, om[k]));

sumKCxiOverI(s, om) =
{
  my(mindeg = vecmin(vector(#om, k, om[k] - 1)));
  my(maxi = N \ mindeg);
  sum(i = 1, maxi, KCsum(subst(s, x, x^i), om) / i);
}

solve_s(om) =
{
  my(s = x + O(x^(N + 2)));
  for(iter = 1, N + 2, s = x * exp(sumKCxiOverI(s, om)));
  s;
}

ZDm(s, m) =
{
  my(ds = divisors(m));
  my(part1 = sum(k = 1, #ds,
       eulerphi(ds[k]) / (2 * m) * subst(s, x, x^ds[k])^(m / ds[k])));
  my(p1 = s, p2 = subst(s, x, x^2));
  my(part2 = if(m % 2 == 1,
       1/2 * p1 * p2^((m - 1) / 2)
     ,
       1/4 * (p1^2 * p2^((m - 2) / 2) + p2^(m / 2))
     ));
  part1 + part2;
}

solve_G(om) =
{
  my(s = solve_s(om));
  my(KC = KCsum(s, om));
  my(E = exp(sumKCxiOverI(s, om)));
  my(SX = x * (E - 1));
  my(SC = (E - 1) - KC);
  my(TS = x * SC);
  my(TSCm = KC * SX);
  my(TCm = sum(k = 1, #om, ZDm(s, om[k])));
  [s, TCm + TS - TSCm];
}

print("Omega = ", Om);
print("Solving s(x) and G(x) natively in PARI/GP up to N=", N, "...");
result = solve_G(Om);
s = result[1];
G = result[2];

print("\nRooted series s(x), first non-zero terms:");
rooted_terms = List();
for(n = 0, N, {
  c = polcoeff(s, n);
  if(c != 0, listput(rooted_terms, [n, c]));
});
print(rooted_terms);

print("\nUnrooted series G(x), first non-zero terms:");
unrooted_terms = List();
for(n = 0, N, {
  c = polcoeff(G, n);
  if(c != 0, listput(unrooted_terms, [n, c]));
});
print(unrooted_terms);

print("\nCross-check against mgonal_cactus_series_omega.py (Python,");
print("fractions.Fraction) and verify_dissymmetry_omega.py (Python,");
print("sympy.Rational) for Omega={5,6}:");
published_rooted = [[1,1],[5,1],[6,1],[9,3],[10,6],[11,4],[13,13],[14,41],[15,49],[16,22]];
published_unrooted = [[5,1],[6,1],[9,1],[10,1],[11,1],[13,3],[14,6],[15,6]];
print("  expected rooted (first 10): ", published_rooted);
print("  expected unrooted (first 8): ", published_unrooted);

match_rooted = 1;
for(k = 1, #published_rooted, {
  my(n = published_rooted[k][1]);
  my(expected = published_rooted[k][2]);
  if(polcoeff(s, n) != expected, match_rooted = 0);
});
match_unrooted = 1;
for(k = 1, #published_unrooted, {
  my(n = published_unrooted[k][1]);
  my(expected = published_unrooted[k][2]);
  if(polcoeff(G, n) != expected, match_unrooted = 0);
});
print("\nRooted match: ", match_rooted == 1);
print("Unrooted match: ", match_unrooted == 1);

quit;
