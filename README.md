# Verification code for "Enumeration and Asymptotic Analysis of Strict Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths"

**→ [See the companion page](docs/index.html)** - the four Table 1 examples, a growth-rate chart, and the eight sequences prepared for OEIS, illustrated.

## Rationale

In [1], I enumerate strict non-plane *m*-gonal cactus graphs for a single fixed cycle length *m* ≥ 5. Its concluding remarks note that the method extends "without difficulty" to a finite mixed set Ω of admissible cycle lengths, but this extension was never carried out there, for any Ω, numerically or analytically.

This repository accompanies the paper that carries it out [2]. For an arbitrary finite Ω, the paper gives an exact characterization of which vertex counts occur (reducing, for two admissible lengths, to the classical Frobenius coin problem), a closed-form critical value when every length in Ω is odd, a proof that this same method is structurally obstructed as soon as an even length is present, and the general asymptotic enumeration law for arbitrary finite Ω, with a closed form for an associated second-order coefficient in the all-odd case.

The preprint of [1] is available on Zenodo (DOI [10.5281/zenodo.21513753](https://zenodo.org/records/21513753)), and its accompanying code on GitHub: [non-plane-mgonal-cacti](https://github.com/frederic-speyser/non-plane-mgonal-cacti). Paper [2], titled "Enumeration and Asymptotic Analysis of Strict Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths," is in preparation for submission to the Journal of Integer Sequences. A working-paper version and the preprint, presenting the same results, will be placed in References.

## Related repositories

- **[non-plane-mgonal-cacti](https://github.com/frederic-speyser/non-plane-mgonal-cacti)** - code and data for [1], the single-length enumeration paper this work generalizes. `mgonal_cactus_series_omega.py` below is a direct generalization of `mgonal_cactus_series.py` from that repository; no other code is shared between the two.
- **[cactus-growth-rate-monotonicity](https://github.com/frederic-speyser/cactus-growth-rate-monotonicity)** - proves the growth-rate monotonicity conjecture of [1] for a single length *m*. Unrelated in content: that paper explicitly does not address the mixed-Ω case treated here.
- **[cactus-split-decomp-omega](https://github.com/frederic-speyser/cactus-split-decomp-omega)** - an earlier, exploratory, non-peer-reviewed numerical investigation of the single case Ω = {5, 6}, predating and motivating paper [2]. Cited in [2] as preliminary work; superseded, for Ω = {5, 6}, by the exact and proved results here.

## Repository layout

```
├── python/               the main solvers, plus 5 independent cross-check scripts
├── pari/                 one more independent cross-check, in PARI/GP rather than Python
├── tests/                automated regression tests (pytest)
├── bfiles/               100-term OEIS b-files for the 8 sequences below
├── docs/                 a companion illustrative page (index.html)
├── CHANGELOG.md
├── LICENSE
└── README.md             (this file / ce fichier)
```

## Main pipeline (`python/`)

| File | What it computes |
|---|---|
| `mgonal_cactus_series_omega.py` | Exact rooted and unrooted enumeration series, indexed by **vertex count**, generalizing `mgonal_cactus_series.py` from [1]: the kernel *K<sub>C</sub>* becomes a sum of one term per size in Ω. Exact rational arithmetic (Python `Fraction`). Used to produce Table 1 of paper [2] (ρ_Ω, τ_Ω, and the exact support characterization of Theorem 3). |
| `mgonal_cactus_series_omega_blocks.py` | Exact rooted and unrooted enumeration series, indexed by **number of blocks** (the convention this OEIS family already uses). Evaluates the block-indexed functional equation directly, rather than tracking vertex count and block count separately - computes 100 terms per sequence in about a minute (up to ~70s for the heaviest Ω tested, Ω={5,7,9}). This is what generated the data in `bfiles/`. |
| `critical_point_solver.py` | High-precision direct solver for the critical pair (ρ_Ω, τ_Ω), by a damped fixed-point iteration with Aitken extrapolation. Implements the closed form of Theorem 5 (Ω all-odd) and the consistency check of Proposition 6 (Ω containing an even length). |
| `reproduce_table1.py` | Driver script: reproduces Table 1 of paper [2] end to end, combining the two solvers above. |

## Independent cross-checks (`python/` and `pari/`)

Each of these re-derives a piece of paper [2]'s results by a genuinely different route - different arithmetic, a different construction, or a different language entirely - so that agreement with the main pipeline is evidence from outside that pipeline, not a restatement of it. The four Ω printed in Table 1 of the paper ({5,6}, {5,7}, {5,7,9}, {5,6,7}) are illustrative examples, not the full extent of what is checked: the files marked *(extended)* below run the same kind of check on eight further Ω that do not appear in Table 1, chosen to cover the same qualitative combinations (both sizes |Ω|=2, |Ω|=3 and |Ω|=4, both the all-odd and mixed-parity regimes, and both one and two even cycle lengths present at once)..

**All 12 Ω are gathered in the table below:**

| Ω | \|Ω\| | parity | source | critical-value route (Thm. 5 / Prop. 6) |
|---|---|---|---|---|
| {5,6} | 2 | mixed | Table 1 | Prop. 6 |
| {5,7} | 2 | odd | Table 1 | Thm. 5 |
| {5,7,9} | 3 | odd | Table 1 | Thm. 5 |
| {5,6,7} | 3 | mixed | Table 1 | Prop. 6 |
| {7,11} | 2 | odd | extended | Thm. 5 |
| {6,7} | 2 | mixed | extended | Prop. 6 |
| {5,8} | 2 | mixed | extended | Prop. 6 |
| {7,9,11} | 3 | odd | extended | Thm. 5 |
| {5,9,13} | 3 | odd | extended | Thm. 5 |
| {5,6,8} | 3 | mixed (one even length) | extended | Prop. 6 |
| {5,7,9,11} | 4 | odd | extended | Thm. 5 |
| {6,8,9} | 3 | mixed (two even lengths) | extended | Prop. 6 |

Every Ω here is checked by at least 2 independent implementations (Theorem 3's support characterization plus the two series solvers); the 8 *extended* rows additionally get the PARI/GP cross-check, and the numeric critical-pair consistency check against Theorem 5 or Proposition 6. This is not, and cannot be, a test of "all" finite Ω - there are infinitely many - but between the sizes, parities, and (for the mixed cases) the number of even lengths present, every qualitative combination the theorems distinguish is now exercised at least once.

| File | Independent route |
|---|---|
| `growth_rate_omega.py` | Estimates ρ_Ω from a plain ratio of consecutive series coefficients, instead of solving the criticality condition analytically. |
| `exhaustive_iso_omega.py` | For the smallest cases (1–2 blocks) of Ω={5,6}, builds the cactus graphs directly with `networkx` and counts isomorphism classes by hand — no functional equation at all. |
| `exhaustive_iso_omega_extended.py` *(extended)* | The same direct-construction check, on Ω={7,11} instead. |
| `split_tree_omega.py` | Checks by brute force that the split-decomposition characterization (Theorem 1) still holds when two *different* cycle sizes from Ω={5,6} meet at the same cut vertex — the genuinely new case paper [2] introduces. |
| `split_tree_omega_extended.py` *(extended)* | The same brute-force check, on Ω={7,11} instead. |
| `verify_dissymmetry_omega.py` | Recomputes the unrooted series from scratch using `sympy.Rational` and hand-written convolutions, instead of `fractions.Fraction`. Runs any Ω via `--omega`; used on Table 1's four Ω. |
| `verify_pari_omega.gp` | A third full implementation, in PARI/GP instead of Python, using PARI's native truncated power-series arithmetic, for Ω={5,6}. |
| `verify_extended_omega.py` *(extended)* | For eight Ω not in Table 1: Theorem 3's support check, agreement between the two independent Python series implementations, and consistency of the numerically located critical pair with Theorem 5 / Proposition 6. |
| `verify_pari_omega_extended.gp` *(extended)* | The same eight Ω through PARI/GP, cross-checked against `verify_extended_omega.py`'s output. |

## Tests (`tests/`)

- `test_regression_known_values.py` - solver output anchored against values already published in [1] (Theorem 2 and Table 3), covering both parities of the kernel.
- `test_cross_consistency.py` - agreement between the exact-series route and the direct numerical solver.
- `test_theorem_a_untabulated_omega.py` - paper [2]'s support-characterization theorem (Theorem 3), verified on Ω not appearing in its published Table 1.
- `test_untabulated_omega_extended.py` - the two fast checks of `verify_extended_omega.py` (Theorem 3, cross-implementation agreement), run automatically for eight Ω not appearing in Table 1.
- `test_random_omega_property.py` - Theorem 3's support characterization, on 30 finite sets Ω drawn at random (fixed seed, |Ω| ∈ {2,3,4}) rather than chosen by hand, to guard against unintentional selection bias in every other Ω tested in this repository.
- `test_boundary_single_length_omega.py` - the |Ω|=1 boundary case: checks that the general block-indexed series, applied to a singleton Ω, reproduces series independently published and verified on the OEIS (A398033, A397250, A397210) - not just the critical pair already anchored by `test_regression_known_values.py`, but the full series.
- `test_vertex_indexed_regression.py` - locks the vertex-indexed engine (`mgonal_cactus_series_omega.py`) against silent regressions.
- `test_oeis_data_regression.py` - locks the eight sequences prepared for OEIS submission (block-indexed) against silent regressions.

## Method note: a direct univariate formulation for the block-indexed series

The block-indexed series can be obtained 2 ways: by tracking vertex count and block count as 2
separate variables and summing out the vertex dimension at the end, or by evaluating the same
functional equation directly in the block-counting variable alone, since for each fixed number of
blocks only finitely many vertex counts contribute; the 2 are the same power series identity,
just reached by different routes. `mgonal_cactus_series_omega_blocks.py` uses the direct route: it
is substantially faster, since it never needs to track the vertex dimension at all, which is what
made the bivariate route expensive at high term counts. Both routes were implemented and
cross-checked against each other, and against the independently published data of the companion
exploratory repository for Ω={5,6}, before the direct route became the one used to generate this
repository's data.

## Requirements

Python 3, with [`mpmath`](https://mpmath.org/) (`pip install mpmath`) and, for two of the cross-checks, [`sympy`](https://www.sympy.org/) and [`networkx`](https://networkx.org/) (`pip install sympy networkx`). [`pytest`](https://pytest.org/) is optional but recommended for the test suite (`pip install pytest`). The PARI/GP cross-check additionally needs a `gp` installation (see [pari.math.u-bordeaux.fr](https://pari.math.u-bordeaux.fr/)).

Exact versions this repository's checks were last run against, for strict reproducibility (any reasonably recent version of each should work; these are not hard pins, just the record of what was actually used):

| Package | Version |
|---|---|
| Python | 3.11.15 |
| mpmath | 1.3.0 |
| sympy | 1.14.0 |
| networkx | 3.6.1 |
| pytest | 9.1.1 |

## Usage

```bash
python3 python/reproduce_table1.py
python3 -m pytest tests/ -v

# independent cross-checks (Table 1's four Omega)
python3 python/growth_rate_omega.py
python3 python/exhaustive_iso_omega.py
python3 python/split_tree_omega.py
python3 python/verify_dissymmetry_omega.py
gp -q pari/verify_pari_omega.gp < /dev/null

# extended independent cross-checks (eight further, untabulated Omega)
python3 python/verify_extended_omega.py   # several minutes; not part of pytest
python3 python/exhaustive_iso_omega_extended.py
python3 python/split_tree_omega_extended.py
gp -q pari/verify_pari_omega_extended.gp < /dev/null
```

## Data availability and OEIS submission

There are 4 mixed sets Ω covered here — {5,6}, {5,7}, {5,7,9}, {5,6,7}. For each one, Table 1 of paper [2] gives a two-number summary of its growth rate (ρ_Ω, τ_Ω, and whether a closed form exists) - 1 row per Ω, 4 rows in total. Separately, for each of those same four Ω, the *full* rooted and unrooted term-by-term counts (1, 2, 13, 125, 1393, ...) are computed and prepared for OEIS - 2sequences per Ω, 8 sequences in total. The four-row summary and the eight full sequences describe the same 4 cases at two different levels of detail, not 8 different cases; the [companion page](docs/index.html) keeps them in two separate tables for that reason, rather than mixing a growth-rate summary with raw term data in one table.

All 8 sequences were computed by `mgonal_cactus_series_omega_blocks.py`, indexed by number of blocks per the convention this OEIS family already uses. 100 verified terms per sequence are provided in `bfiles/`. See [oeis.org/search?q=speyser](https://oeis.org/search?q=speyser) for current submission status and A-numbers once assigned.

What makes them new:

- **The mixture itself.** Every non-plane-cactus sequence in this family submitted so far - A332648/A332649 (general arrays, by Andrew Howroyd) and the single-length columns I submitted separately for paper [1] - fixes 1 cycle length *m* throughout. These 8 are, as far as I can tell, the first tabulated counts for a *mix* of at least two different cycle lengths at once, which is exactly what paper [2] adds and the earlier papers did not attempt.
- **A closed form where there wasn't one before.** For the all-odd cases (Ω = {5,7} and Ω = {5,7,9}), the growth rate is given by an exact closed form (Theorem 5 of [2]), not just a numerically located constant.
- **An exact description of where the sequence is nonzero.** For the vertex-indexed version of each sequence, the full set of nonzero positions is characterized (Theorem 3 of [2]) via the numerical semigroup generated by {m − 1 : m ∈ Ω} - reducing, for two generators, to the classical Frobenius coin problem.


## References

[1] Speyser, F. G. *Enumeration and Asymptotic Analysis of Strict Non-Plane m-Gonal Cactus Graphs via Split-Decomposition.* Submitted to the Electronic Journal of Combinatorics, 2026. Preprint: DOI [10.5281/zenodo.21513753](https://zenodo.org/records/21513753).

[2] Speyser, F. G. *Enumeration and Asymptotic Analysis of Strict Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths.* Working paper, in preparation for submission to the Journal of Integer Sequences, 2026.

## Citation

If you use this code, please cite the papers above. A citable archive of this repository is available via Zenodo: DOI: xxxxx.

## Author

Frédéric G. Speyser - Independent Researcher, France
ORCID: [0000-0002-1767-5325](https://orcid.org/0000-0002-1767-5325)

## License

MIT (see `LICENSE`), for consistency with the other repositories.
