"""
test_oeis_sequences_regression.py

Locks mgonal_cactus_series_omega.py against silent regressions on the
eight sequences already prepared for submission to the OEIS (see
oeis_submissions.md, alongside this repository). Once these sequences
are actually submitted and assigned A-numbers, any future change to this
project's code that altered their values would otherwise go unnoticed
until caught by an OEIS editor or a reader -- this test catches it here
instead, before it ever leaves the repository.

Data below reproduced exactly from oeis_submissions.md (with explicit
zeros, matching the convention adopted there), truncated to the first 35
terms of each sequence (a(1) through a(35)).

Reference: F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths" (in preparation 
for submission to the Journal of Integer Sequences, 2026); oeis_submissions.md.

Author: Frederic G. Speyser
Run: python3 -m pytest test_oeis_sequences_regression.py -v
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import mgonal_cactus_series_omega as om

EXPECTED = {
    (5, 6): {
        'rooted': [1,0,0,0,1,1,0,0,3,6,4,0,13,41,49,22,62,278,498,415,473,1920,
                   4600,5693,5547,14359,40326,66324,74199,126743,349403,703367,
                   953136,1345423,3134905],
        'unrooted': [0,0,0,0,1,1,0,0,1,1,1,0,3,6,6,4,8,25,42,32,44,140,302,357,
                     353,848,2192,3391,3759,6300,16348,31201,41008,57286,128671],
    },
    (5, 7): {
        'rooted': [1,0,0,0,1,0,1,0,3,0,6,0,17,0,44,0,116,0,331,0,917,0,2673,0,
                   7709,0,22742,0,67332,0,201125,0,604950,0,1828308],
        'unrooted': [0,0,0,0,1,0,1,0,1,0,1,0,4,0,6,0,14,0,30,0,76,0,185,0,485,
                     0,1273,0,3471,0,9524,0,26665,0,75244],
    },
    (5, 7, 9): {
        'rooted': [1,0,0,0,1,0,1,0,4,0,6,0,24,0,52,0,178,0,460,0,1503,0,4287,0,
                   13761,0,41378,0,132433,0,411116,0,1319591,0,4179208],
        'unrooted': [0,0,0,0,1,0,1,0,2,0,1,0,5,0,7,0,22,0,40,0,123,0,289,0,858,
                     0,2261,0,6737,0,19109,0,57390,0,169218],
    },
    (5, 6, 7): {
        'rooted': [1,0,0,0,1,1,1,0,3,6,10,7,17,41,93,120,178,342,829,1469,2332,
                   3828,8121,16358,29543,49319,92160,184231,358482,642696,
                   1159976,2211938,4348394,8227717,15136332],
        'unrooted': [0,0,0,0,1,1,1,0,1,1,2,1,4,6,12,13,21,32,72,109,177,268,
                     545,1000,1744,2762,5007,9437,17647,30326,53213,97707,
                     185385,338123,604179],
    },
}


def test_oeis_sequences_unchanged():
    for omega, expected in EXPECTED.items():
        om.N = 45
        s, G = om.solve_G(omega)
        rooted = [c.numerator for n, c in enumerate(s) if n >= 1][:35]
        unrooted = [c.numerator for n, c in enumerate(G) if n >= 1][:35]
        assert rooted == expected['rooted'], f"Omega={omega}: rooted sequence changed"
        assert unrooted == expected['unrooted'], f"Omega={omega}: unrooted sequence changed"


if __name__ == "__main__":
    test_oeis_sequences_unchanged()
    print(f"All {len(EXPECTED) * 2} OEIS-submission sequences reproduced exactly.")
