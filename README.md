# Verification code for "Enumeration and Asymptotic Analysis of Strict Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths"


## Rationale

In [1], I enumerate strict non-plane *m*-gonal cactus graphs for a single fixed cycle length *m* ≥ 5. Its concluding remarks note that the method extends "without difficulty" to a finite mixed set Ω of admissible cycle lengths, but this extension was never carried out there, for any Ω, numerically or analytically.

This repository accompanies the paper that carries it out. For an arbitrary finite Ω, the paper gives an exact characterization of which vertex counts occur (reducing, for two admissible lengths, to the classical Frobenius coin problem), a closed-form critical value when every length in Ω is odd, a proof that this same method is structurally obstructed as soon as an even length is present, and the general asymptotic enumeration law for arbitrary finite Ω, with a closed form for an associated second-order coefficient in the all-odd case.

The preprint of [1] is available on Zenodo (DOI [10.5281/zenodo.21513753](https://doi.org/10.5281/zenodo.21513753)), and its accompanying code on GitHub: [non-plane-mgonal-cacti](https://github.com/frederic-speyser/non-plane-mgonal-cacti). This paper — the one this repository provides verification code for — has been submitted to the *Journal of Integer Sequences*.

## Related repositories

- **[non-plane-mgonal-cacti](https://github.com/frederic-speyser/non-plane-mgonal-cacti)** - code and data for [1], the single-length enumeration paper this work generalizes. `mgonal_cactus_series_omega.py` below is a direct generalization of `mgonal_cactus_series.py` from that repository; no other code is shared between the two.
- **[cactus-growth-rate-monotonicity](https://github.com/frederic-speyser/cactus-growth-rate-monotonicity)** - proves the growth-rate monotonicity conjecture of [1] for a single length *m*. Unrelated in content: that paper explicitly does not address the mixed-Ω case treated here.
- **[cactus-split-decomp-omega](https://github.com/frederic-speyser/cactus-split-decomp-omega)** - an earlier, exploratory, non-peer-reviewed numerical investigation of the single case Ω = {5, 6}, predating and motivating the present paper. Cited in [1's continuation] as preliminary work; superseded, for Ω = {5, 6}, by the exact and proved results here.

## What this repository contains

- **`mgonal_cactus_series_omega.py`** - exact rooted and unrooted enumeration series for strict cactus graphs admitting a finite set Ω of cycle lengths, generalizing `mgonal_cactus_series.py` from [1]: the kernel *K<sub>C</sub>* becomes a sum of one USEQ term per size in Ω. Exact rational arithmetic (Python `Fraction`). This is the script that generated the enumerative data tabulated in the paper and submitted to the OEIS.
- **`critical_point_solver.py`** - high-precision direct solver for the critical pair (ρ_Ω, τ_Ω), by a damped fixed-point iteration with Aitken extrapolation. Implements the closed form of Theorem 5 (Ω all-odd) and the consistency check of Proposition 6 (Ω containing an even length), and the closed-form second-order coefficient of Theorem 8(b).
- **`reproduce_table1.py`** — driver script: reproduces the paper's Table 1 end to end, combining the two modules above.
- **`tests/`** — four automated tests (pytest-compatible):
  - `test_regression_known_values.py` — solver output anchored against values already published in [1] (Theorem 2 and Table 3), covering both parities of the kernel.
  - `test_cross_consistency.py` - agreement between the exact-series route and the direct numerical solver.
  - `test_theorem_a_untabulated_omega.py` - the paper's support-characterization theorem, verified on Ω not appearing in the published Table 1.
  - `test_oeis_sequences_regression.py` - locks the eight sequences prepared for OEIS submission against silent future regressions.

## Requirements

Python 3, with `mpmath` (`pip install mpmath`). `pytest` is optional but recommended for the test suite (`pip install pytest`).

## Usage

```bash
python3 reproduce_table1.py
python3 -m pytest tests/ -v
```

## Data availability

Eight integer sequences (rooted and unrooted, for Ω = {5,6}, {5,7}, {5,7,9}, {5,6,7}) have been computed by `mgonal_cactus_series_omega.py` and prepared for submission to the OEIS; see the paper for the full data and the current submission status.

## References

[1] Speyser, F. G. *Enumeration and Asymptotic Analysis of Strict Non-Plane m-Gonal Cactus Graphs via Split-Decomposition.* Submitted to the Electronic Journal of Combinatorics, 2026. Preprint: DOI [10.5281/zenodo.21513753](https://doi.org/10.5281/zenodo.21513753).

## Citation

If you use this code, please cite the paper above. A citable archive of this repository will be added here via Zenodo once deposited.

## Author

Frédéric G. Speyser - Independent Researcher, France - ORCID: [0000-0002-1767-5325](https://orcid.org/0000-0002-1767-5325)

## License

MIT (see `LICENSE`), for consistency with the other repositories.
