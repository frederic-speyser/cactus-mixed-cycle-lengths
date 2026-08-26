"""
test_boundary_single_length_omega.py

A boundary-case regression test: |Omega| = 1. Every other check in this
repository exercises the genuinely new case paper [2] introduces -- a
*mixed* set of at least two cycle lengths. None of them directly confirms
that the general machinery, when handed a singleton Omega = {m}, reduces
exactly to the single-length construction of the companion paper [1] --
the case the general code is explicitly built to specialize to (see the
module docstrings of mgonal_cactus_series_omega.py and
mgonal_cactus_series_omega_blocks.py, and Theorem 3/5's own text, which
state this reduction but do not, by themselves, test it).

test_regression_known_values.py already anchors the *critical pair*
(rho_m, tau_m) for two singletons (m=5, m=6) against values published in
[1]. This file complements that by anchoring the full block-indexed
*series* -- not just its critical point -- for singleton Omega against
data independently published and verified on the OEIS, checked directly
against the live OEIS pages on 2026-08-26 (not against any value computed
internally by this project, so a bug shared by this repository's own
series solver and its own cross-checks could not silently pass here):

  - Omega = (5,): rooted terms match A398033 ("Number of rooted unlabeled
    5-gonal cacti having n blocks"), offset 0.
  - Omega = (5,): unrooted terms match A397250 ("Number of unlabeled
    5-gonal cacti having n blocks"), offset 1.
  - Omega = (7,): rooted terms match A397210 ("Number of rooted unlabeled
    7-gonal cacti having n blocks"), offset 0.

(A398575, the k=7 unrooted entry, is allocated on OEIS but not yet
populated with data as of this writing, so it cannot serve as a ground
truth here; it is omitted rather than compared against a placeholder.)

Reference:
  [1] F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
      Non-Plane m-Gonal Cactus Graphs via Split-Decomposition."
  [2] F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
      Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths" (in
      preparation for submission to the Journal of Integer Sequences,
      2026).

Author: Frederic G. Speyser
Run: python3 -m pytest tests/test_boundary_single_length_omega.py -v
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))

from mgonal_cactus_series_omega_blocks import compute

# Hardcoded from the live OEIS pages, checked 2026-08-26 -- not derived
# from this repository's own code.
A398033_ROOTED_M5 = [
    1, 1, 3, 13, 62, 333, 1894, 11258, 68990, 432964, 2767569, 17957046,
    117951892, 782810692, 5241104161, 35357170219,
]
A397250_UNROOTED_M5 = [
    1, 1, 3, 8, 31, 132, 636, 3280, 17958, 101877, 595796, 3564222,
    21731034, 134586653, 844796689,
]
A397210_ROOTED_M7 = [
    1, 1, 4, 25, 176, 1397, 11757, 103376, 937179, 8699140, 82250503,
]


def test_omega5_rooted_matches_A398033():
    rooted, _ = compute((5,), terms=len(A398033_ROOTED_M5) - 1)
    assert rooted == A398033_ROOTED_M5


def test_omega5_unrooted_matches_A397250():
    _, unrooted = compute((5,), terms=len(A397250_UNROOTED_M5))
    assert unrooted[:len(A397250_UNROOTED_M5)] == A397250_UNROOTED_M5


def test_omega7_rooted_matches_A397210():
    rooted, _ = compute((7,), terms=len(A397210_ROOTED_M7) - 1)
    assert rooted == A397210_ROOTED_M7


if __name__ == "__main__":
    test_omega5_rooted_matches_A398033()
    test_omega5_unrooted_matches_A397250()
    test_omega7_rooted_matches_A397210()
    print("Singleton Omega reduction verified against A398033, A397250, and A397210.")
