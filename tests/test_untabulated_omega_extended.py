"""
test_untabulated_omega_extended.py

Automated regression test wrapping the two FAST checks of
verify_extended_omega.py (Theorem 3's support characterization, and
agreement between the two independent Python series implementations)
for all six untabulated Omega, so they run automatically under
`pytest tests/ -v` like the rest of the suite.

The third check in verify_extended_omega.py -- numeric location of the
critical pair and its consistency with Theorem 5 / Proposition 6 -- is
deliberately NOT included here: it takes on the order of a minute or two
per Omega at the precision needed for a meaningful check (see that
script's own comments), which would make the default test suite far
slower without adding a different kind of evidence than what is already
covered by test_cross_consistency.py (for the four tabulated Omega) and
reproduce_table1.py. Run verify_extended_omega.py directly for that
slower, more expensive check.

Reference: F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths" (in
preparation for submission to the Journal of Integer Sequences, 2026).

Author: Frederic G. Speyser
Run: python3 -m pytest test_untabulated_omega_extended.py -v
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))

from verify_extended_omega import (
    UNTABULATED_OMEGA, check_theorem3, check_cross_implementation,
)


def test_theorem3_on_all_untabulated_omega():
    for omega in UNTABULATED_OMEGA:
        ok, actual, predicted = check_theorem3(omega)
        assert ok, (
            f"Omega={omega}: Theorem 3 predicts support {sorted(predicted)}, "
            f"actual series support is {sorted(actual)}"
        )


def test_cross_implementation_on_all_untabulated_omega():
    for omega in UNTABULATED_OMEGA:
        ok, s_main, G_main, s_vd, G_vd = check_cross_implementation(omega)
        assert ok, (
            f"Omega={omega}: main solver and sympy solver disagree.\n"
            f"  s(x) [main]:  {s_main}\n  s(x) [sympy]: {s_vd}\n"
            f"  G(x) [main]:  {G_main}\n  G(x) [sympy]: {G_vd}"
        )


if __name__ == "__main__":
    test_theorem3_on_all_untabulated_omega()
    test_cross_implementation_on_all_untabulated_omega()
    print(f"All fast checks passed on {len(UNTABULATED_OMEGA)} untabulated Omega.")
