/*
verify_pari_omega_extended.gp

The PARI/GP analogue of verify_extended_omega.py: runs the same batch of
eight finite sets Omega that do NOT appear in Table 1 of the paper through
PARI's own native truncated power-series arithmetic -- a third,
completely independent implementation (different language, different
runtime, different internal representation of O(x^N) series) alongside
the two Python routes already cross-checked against each other in
verify_extended_omega.py (mgonal_cactus_series_omega.py, using
fractions.Fraction, and verify_dissymmetry_omega.py, using
sympy.Rational).

The functional equation implemented here is identical to
verify_pari_omega.gp (which only ever ran Omega={5,6}); this script just
loops the same machinery over eight untabulated Omega instead (the last
two, {5,7,9,11} and {6,8,9}, added to cover |Omega|=4 and a mixed-parity
Omega with two even lengths simultaneously, neither of which arose in the
first six), and prints the leading terms of both the rooted series s(x)
and the unrooted series G(x) for each, together with a PASS/FAIL
comparison against the terms independently computed by
verify_extended_omega.py (hardcoded below, from that script's own printed
output -- see CHANGELOG.md for the run this was taken from).

Development notes (see verify_pari_omega.gp for the original set):
  - a function body already wrapped in {...} cannot contain a further
    {...}-grouped multi-statement loop; every helper below is written as
    a single expression;
  - the two-argument arrow form vector(n, i -> expr) is NOT accepted by
    this parser (it collides with the three-argument vector(n, i, expr)
    form) -- every vector() call below uses the three-argument form;
    select(pred -> ..., v), by contrast, does accept the arrow form and
    is used that way below.

Reference: F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths" (in
preparation for submission to the Journal of Integer Sequences, 2026).

Author: Frederic G. Speyser
Run: gp -q verify_pari_omega_extended.gp < /dev/null
*/

N = 30;
default(seriesprecision, N + 2);

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

all_nonzero_terms(ser) =
{
  my(degs = select(n -> polcoeff(ser, n) != 0, vector(N + 1, n, n - 1)));
  vector(#degs, k, [degs[k], polcoeff(ser, degs[k])]);
}

first_k(v, k) = vector(min(k, #v), i, v[i]);

OMEGA_LIST = [[7,11], [6,7], [5,8], [7,9,11], [5,9,13], [5,6,8], [5,7,9,11], [6,8,9]];
EXPECTED_ROOTED = [[[1,1],[7,1],[11,1],[13,4],[17,9],[19,25],[21,6],[23,99]], [[1,1],[6,1],[7,1],[11,4],[12,7],[13,4],[16,22],[17,62]], [[1,1],[5,1],[8,1],[9,3],[12,7],[13,13],[15,5],[16,53]], [[1,1],[7,1],[9,1],[11,1],[13,4],[15,8],[17,14],[19,35]], [[1,1],[5,1],[9,4],[13,21],[17,133],[21,937],[25,7050],[29,55318]], [[1,1],[5,1],[6,1],[8,1],[9,3],[10,6],[11,4],[12,7]], [[1,1],[5,1],[7,1],[9,4],[11,7],[13,24],[15,60],[17,187]], [[1,1],[6,1],[8,1],[9,1],[11,4],[13,8],[14,8],[15,5]]];
EXPECTED_UNROOTED = [[[7,1],[11,1],[13,1],[17,1],[19,4],[21,1],[23,9],[25,14]], [[6,1],[7,1],[11,1],[12,1],[13,1],[16,4],[17,7],[18,7]], [[5,1],[8,1],[9,1],[12,1],[13,3],[15,1],[16,7],[17,8]], [[7,1],[9,1],[11,1],[13,1],[15,1],[17,2],[19,5],[21,9]], [[5,1],[9,2],[13,5],[17,17],[21,80],[25,454],[29,2878]], [[5,1],[6,1],[8,1],[9,1],[10,1],[11,1],[12,1],[13,4]], [[5,1],[7,1],[9,2],[11,2],[13,5],[15,8],[17,23],[19,49]], [[6,1],[8,1],[9,1],[11,1],[13,1],[14,1],[15,1],[16,5]]];

all_ok = 1;
for(idx = 1, #OMEGA_LIST, {
  my(om = OMEGA_LIST[idx]);
  print("Omega = ", om);
  my(result = solve_G(om));
  my(s = result[1]);
  my(G = result[2]);
  my(got_rooted = first_k(all_nonzero_terms(s), #EXPECTED_ROOTED[idx]));
  my(got_unrooted = first_k(all_nonzero_terms(G), #EXPECTED_UNROOTED[idx]));
  print("  rooted   (PARI/GP):   ", got_rooted);
  print("  rooted   (expected):  ", EXPECTED_ROOTED[idx]);
  print("  unrooted (PARI/GP):   ", got_unrooted);
  print("  unrooted (expected):  ", EXPECTED_UNROOTED[idx]);
  my(ok = (got_rooted == EXPECTED_ROOTED[idx]) && (got_unrooted == EXPECTED_UNROOTED[idx]));
  print("  match: ", ok);
  if(!ok, all_ok = 0);
  print();
});

print("ALL OMEGA MATCH: ", all_ok);
quit;
