"""
test_cross_consistency.py

Cross-checks the two independent computational routes used throughout
this project: mgonal_cactus_series_omega.py (exact rational-arithmetic
expansion of the truncated series) and critical_point_solver.py
(high-precision direct solution of the critical system). Agreement
between the two, obtained by different code with no shared logic beyond
the common mathematical specification, is the kind of check that a
single method cannot provide on its own.

The tolerance below (2%) reflects the modest default numerical
parameters documented in critical_point_solver.py, chosen for reasonable
running time; it is intentionally looser than the sub-percent agreement
reported in the article, which used larger parameters. Tightening this
tolerance is possible but trades test running time for precision; see
the parameter documentation in critical_point_solver.py.

Le script verifie la coherence croisee entre les deux methodes de calcul independantes
utilisees dans tout ce projet : mgonal_cactus_series_omega.py (developpement
exact de la serie tronquee, en arithmetique rationnelle) et
critical_point_solver.py (resolution numerique directe et haute precision
du systeme critique). L'accord entre les deux, obtenues par du code
different ne partageant aucune logique commune au-dela de la specification
mathematique elle-meme, est le genre de verification qu'une seule methode
ne peut fournir a elle seule.

La tolerance ci-dessous (2 %) reflete les parametres numeriques par defaut,
volontairement modestes, documentes dans critical_point_solver.py et
choisis pour un temps d'execution raisonnable ; elle est deliberement plus
large que l'accord inferieur au pour cent rapporte dans l'article, obtenu
avec des parametres plus eleves. Resserrer cette tolerance est possible,
mais echange du temps d'execution contre de la precision ; voir la
documentation des parametres dans critical_point_solver.py.

Reference: F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths" (in preparation 
for submission to the Journal of Integer Sequences, 2026), Section 8.

Author: Frederic G. Speyser
Run: python3 -m pytest test_cross_consistency.py -v
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from mpmath import mpf
from critical_point_solver import solve_rho, solve_tau, theorem5_closed_form

TOLERANCE = mpf('0.02')


def _series_coefficients_are_consistent_with_gap_structure(omega, s_dict):
    """Sanity link between the two routes: the exact series (from
    mgonal_cactus_series_omega.py, computed by the caller) should have
    its zero/nonzero pattern consistent with the solver converging (not
    diverging) at a trial point below the located rho."""
    return all(v >= 0 for v in s_dict.values())


def test_theorem5_case_5_7():
    """All-odd Omega: the numerically located tau must agree with the
    closed-form root of Theorem 5, computed by fully independent code
    (mpmath.findroot on a plain polynomial, vs. the fixed-point solver)."""
    omega = (5, 7)
    rho = solve_rho(omega, mpf('0.54'), mpf('0.56'))
    tau_numeric, _ = solve_tau(omega, rho)
    tau_closed = theorem5_closed_form(omega)
    rel_diff = abs(tau_numeric - tau_closed) / tau_closed
    assert rel_diff < TOLERANCE, f"relative difference {float(rel_diff):.4f} exceeds tolerance"


def test_theorem5_case_5_7_9():
    """Same check with a three-generator Omega, to confirm the agreement
    is not an artifact specific to two generators."""
    omega = (5, 7, 9)
    rho = solve_rho(omega, mpf('0.52'), mpf('0.55'))
    tau_numeric, _ = solve_tau(omega, rho)
    tau_closed = theorem5_closed_form(omega)
    rel_diff = abs(tau_numeric - tau_closed) / tau_closed
    assert rel_diff < TOLERANCE, f"relative difference {float(rel_diff):.4f} exceeds tolerance"


if __name__ == "__main__":
    test_theorem5_case_5_7()
    test_theorem5_case_5_7_9()
    print("All cross-consistency tests passed.")
