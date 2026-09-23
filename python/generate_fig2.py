"""
Generates the asymptotic-amplitude figure (rho_Omega-normalized
coefficients s_n compared against the amplitude constant C_Omega) for
the four Omega families listed in Table 1.

For each Omega, plots
    R_n = s_n * rho_Omega^n * n^{3/2}   (n in the admissible class 1+d*k)
together with a horizontal reference line at C_Omega.

rho_Omega and C_Omega plotted here are the values already established by
critical_point_solver.py and amplitude_check_omega.py (REFERENCE_VALUES
below). solve_critical_pair is called to (a) obtain the coefficients s_n
needed for the plot and (b) verify, by explicit assertion, that it
reproduces those already-established values.

Run: python3 generate_fig2.py

Authors: Frederic G. Speyser & Kseniya Vyatkina
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina. All rights reserved.
"""
import mpmath
from mpmath import mp, mpf

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from critical_point_solver import compute_series, solve_critical_pair
from amplitude_check_omega import check as amplitude_check, _d

mp.dps = 30
NMAX = 305

# Starting points for Newton's method only (approximate).
INITIAL_GUESS = {
    (5, 6): ('0.531', '0.712'), (5, 7): ('0.550', '0.727'),
    (5, 7, 9): ('0.535', '0.692'), (5, 6, 7): ('0.506', '0.665'),
}

# Reference values already established by critical_point_solver.py and
# amplitude_check_omega.py; these are the values actually plotted.
REFERENCE_VALUES = {
    (5, 6):     dict(rho='0.530923625274025', C='0.156351357866176'),
    (5, 7):     dict(rho='0.549884575826361', C='0.302088001203884'),
    (5, 7, 9):  dict(rho='0.535218127227685', C='0.269362679524554'),
    (5, 6, 7):  dict(rho='0.506181465494348', C='0.138735005523206'),
}
CONSISTENCY_TOL = mpf('1e-13')

OMEGAS = [(5, 6), (5, 7), (5, 7, 9), (5, 6, 7)]

# Plot styling (kept separate from the numerical computation above).
COLOR_CURVE = '#1f4e79'
COLOR_REFERENCE = '#c0392b'


def dense_Rn(omega, nmin=8):
    d = _d(omega)
    rho, tau, s, N, r1, r2 = solve_critical_pair(omega, *INITIAL_GUESS[omega], N=NMAX)

    ref = REFERENCE_VALUES[omega]
    assert abs(rho - mpf(ref['rho'])) < CONSISTENCY_TOL, \
        f"recomputed rho does not match the reference value for {omega}"
    rho_plot = mpf(ref['rho'])

    ns, Rn = [], []
    k = 0
    while True:
        n = 1 + d * k
        if n > N:
            break
        if n >= nmin and s[n] != 0:
            # The support theorem (aperiodicity of the class 1+dk) is an
            # asymptotic statement: it guarantees s_n != 0 for n large
            # enough, not for every n in the class from the start. The
            # s[n] != 0 filter is therefore necessary for small n; no
            # such zero coefficient occurs beyond n ~ 30 for any of the
            # four families here.
            R = float(s[n] * rho_plot ** n * mpf(n) ** mpf('1.5'))
            ns.append(n)
            Rn.append(R)
        k += 1
    return ns, Rn, float(rho_plot)


if __name__ == "__main__":
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
    for ax, omega in zip(axes.flat, OMEGAS):
        ns, Rn, rho = dense_Rn(omega)
        amp = amplitude_check(omega)
        ref = REFERENCE_VALUES[omega]
        assert abs(amp['C_dir'] - mpf(ref['C'])) < CONSISTENCY_TOL, \
            f"recomputed C_Omega does not match the reference value for {omega}"
        C = float(mpf(ref['C']))

        ax.plot(ns, Rn, '-', color=COLOR_CURVE, linewidth=1.2,
                label=r'$s_n\,\rho_\Omega^n\,n^{3/2}$')
        ax.axhline(C, color=COLOR_REFERENCE, linestyle='--', linewidth=1.2,
                    label=fr'$C_\Omega={C:.6f}$')
        ax.set_title(r'$\Omega=\{' + ','.join(map(str, omega)) + r'\}$')
        ax.set_xlabel(r'$n$')
        ax.set_ylabel(r'$s_n\rho_\Omega^n n^{3/2}$')
        ax.legend(fontsize=8, loc='upper right')
        ax.grid(alpha=0.25)
    fig.suptitle(r'Convergence of $s_n\rho_\Omega^n n^{3/2}$ to $C_\Omega$', fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig('fig2_asymp_corrected.pdf')
    fig.savefig('fig2_asymp_corrected.png', dpi=200)
    fig.savefig('fig2_asymp_corrected.svg')
    print("saved fig2_asymp_corrected.{pdf,png,svg}")
