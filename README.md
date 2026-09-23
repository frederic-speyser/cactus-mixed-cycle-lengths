# Reproducibility and Verification Code for "Enumeration and Asymptotic Analysis of Strict Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths"

[![Tests](https://github.com/frederic-speyser/cactus-mixed-cycle-lengths/actions/workflows/tests.yml/badge.svg)](https://github.com/frederic-speyser/cactus-mixed-cycle-lengths/actions/workflows/tests.yml)

**→ [See the companion page](https://frederic-speyser.github.io/cactus-mixed-cycle-lengths/)** for an illustrated overview: the four Table 1 examples, a growth-rate chart, and the eight sequences.

This repository reproduces all entries of Table 1 of the accompanying paper, and independently cross-checks several of its main computational consequences — the support characterization, the critical-value relations, and the asymptotic amplitude constant — using separate implementations and verification routes.

## Rationale

In *Enumeration and Asymptotic Analysis of Strict Non-Plane m-Gonal Cactus Graphs via Split-Decomposition*, 2026 [1] (DOI [10.5281/zenodo.21513753](https://zenodo.org/records/21513753), and its accompanying code on GitHub: [non-plane-mgonal-cacti](https://github.com/frederic-speyser/non-plane-mgonal-cacti)), F. G. Speyser enumerates strict non-plane *m*-gonal cactus graphs for a single fixed cycle length *m* ≥ 5. Its concluding remarks note that the method extends "without difficulty" to a finite mixed set Ω of admissible cycle lengths, but this extension was never carried out there, for any Ω, numerically or analytically.

This repository accompanies the paper: Speyser, F. G. & Vyatkina, K. *Enumeration and Asymptotic Analysis of Strict Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths*, 2026 [2], DOI: [10.5281/zenodo.22875493](https://doi.org/10.5281/zenodo.22875492). For a finite, nonempty set Ω of cycle lengths ≥ 5, the paper gives an exact characterization of which vertex counts occur (reducing, for two admissible lengths, to the classical Frobenius coin problem), a closed-form critical value when every length in Ω is odd, a proof that the elimination argument giving this closed form does not extend to sets containing an even length, and the general asymptotic form of the counting sequence for arbitrary finite Ω, including an explicit second-derivative contribution entering the leading asymptotic amplitude, in closed form in the all-odd case.

## Related repositories

- **[non-plane-mgonal-cacti](https://github.com/frederic-speyser/non-plane-mgonal-cacti)** - code and data for [1], the single-length enumeration paper this work generalizes. `mgonal_cactus_series_omega.py` below is a direct generalization of `mgonal_cactus_series.py` from that repository; no other code is shared between the two.
- **[cactus-growth-rate-monotonicity](https://github.com/frederic-speyser/cactus-growth-rate-monotonicity)** - proves the growth-rate monotonicity conjecture of [1] for a single length *m*. Unrelated in content: that paper explicitly does not address the mixed-Ω case treated here.
- **[cactus-split-decomp-omega](https://github.com/frederic-speyser/cactus-split-decomp-omega)** - an earlier, exploratory, non-peer-reviewed numerical investigation of the single case Ω = {5, 6}, predating and motivating paper [2]. Cited in [2] as preliminary work; superseded, for Ω = {5, 6}, by the exact and proved results here.

## Quick start

```bash
git clone https://github.com/frederic-speyser/cactus-mixed-cycle-lengths
cd cactus-mixed-cycle-lengths
pip install mpmath pytest
python3 python/reproduce_table1.py
python3 -m pytest tests/ -v
```

This reproduces Table 1 of the paper and runs the full regression suite (14 tests). The independent cross-checks and the PARI/GP scripts need the additional packages listed under Requirements.

## Repository layout

```
├── python/               the main solvers, ten structural/enumeration cross-check scripts,
│                         two asymptotic-amplitude cross-checks, and a figure-generation script
├── pari/                 two more independent cross-checks, in PARI/GP rather than Python
├── tests/                automated regression tests (pytest)
├── bfiles/               100-term OEIS b-files for the eight sequences below
├── docs/                 a companion illustrative page (index.html), for GitHub Pages
├── CHANGELOG.md
├── LICENSE
└── README.md             (this file)
```

## Main pipeline (`python/`)

| File | What it computes |
|---|---|
| `mgonal_cactus_series_omega.py` | Exact rooted and unrooted enumeration series, indexed by **vertex count**, generalizing `mgonal_cactus_series.py` from [1]: the kernel *K<sub>C</sub>* becomes a sum of one term per size in Ω. Exact rational arithmetic (Python `Fraction`). Used to produce Table 1 of paper [2] (ρ_Ω, τ_Ω, and the exact support characterization of Theorem 3). |
| `mgonal_cactus_series_omega_blocks.py` | Exact rooted and unrooted enumeration series, indexed by **number of blocks** (the convention this OEIS family already uses). Evaluates the block-indexed functional equation directly, rather than tracking vertex count and block count separately — computes 100 terms per sequence in about a minute (up to ~70s for the heaviest Ω tested, Ω={5,7,9}). This is what generated the data in `bfiles/` and the OEIS submission drafts. See the *Method note* below for why a direct, univariate route was worth building. |
| `critical_point_solver.py` | High-precision direct solver for the critical pair (ρ_Ω, τ_Ω), by a damped fixed-point iteration with Aitken extrapolation. Implements the closed form of Theorem 10 (Ω all-odd) and the consistency check of Proposition 11 (Ω containing an even length). |
| `reproduce_table1.py` | Driver script: reproduces Table 1 of paper [2] end to end, combining the two solvers above. |
| `generate_fig2.py` | Regenerates the companion page's asymptotic-amplitude figure from the coefficients and constants established by the pipeline above. |

### Method note: a direct univariate formulation for the block-indexed series

The block-indexed series can be obtained two ways: by tracking vertex count and block count as two separate variables and summing out the vertex dimension at the end, or by evaluating the same functional equation directly in the block-counting variable alone, since for each fixed number of blocks only finitely many vertex counts contribute; the two constructions yield the same block-indexed power series, but reach it through different intermediate representations. `mgonal_cactus_series_omega_blocks.py` uses the direct route: it is substantially faster, since it never needs to track the vertex dimension at all, which is what made the bivariate route expensive at high term counts. Both routes were implemented and cross-checked against each other, and against the data previously released in the exploratory repository for Ω={5,6}, before the direct route became the one used to generate this repository's data; see `CHANGELOG.md` for the record of that comparison and for the full history of this repository's reorganization.

## Independent cross-checks (`python/` and `pari/`)

Each of these checks a computational consequence of paper [2]'s results by a separate route — different arithmetic, a different construction, or a different language entirely — so that agreement with the main pipeline is evidence from outside that pipeline, not a restatement of it.

### Coverage across twelve Ω

The four Ω printed in Table 1 of the paper ({5,6}, {5,7}, {5,7,9}, {5,6,7}) are illustrative examples, not the full extent of what is checked: the files marked *(extended)* below run the same kind of check on eight further Ω that do not appear in Table 1, chosen to cover the same qualitative combinations (both sizes |Ω|=2, |Ω|=3 and |Ω|=4, both the all-odd and mixed-parity regimes, and both one and two even cycle lengths present at once). All twelve are gathered here, so the full coverage can be checked at a glance:

| Ω | \|Ω\| | parity | source | critical-value route (Thm. 10 / Prop. 11) |
|---|---|---|---|---|
| {5,6} | 2 | mixed | Table 1 | Prop. 11 |
| {5,7} | 2 | odd | Table 1 | Thm. 10 |
| {5,7,9} | 3 | odd | Table 1 | Thm. 10 |
| {5,6,7} | 3 | mixed | Table 1 | Prop. 11 |
| {7,11} | 2 | odd | extended | Thm. 10 |
| {6,7} | 2 | mixed | extended | Prop. 11 |
| {5,8} | 2 | mixed | extended | Prop. 11 |
| {7,9,11} | 3 | odd | extended | Thm. 10 |
| {5,9,13} | 3 | odd | extended | Thm. 10 |
| {5,7,8} | 3 | mixed (one even length) | extended | Prop. 11 |
| {5,7,9,11} | 4 | odd | extended | Thm. 10 |
| {6,8,9} | 3 | mixed (two even lengths) | extended | Prop. 11 |

Every Ω here is checked by at least two independent computational routes, together with an independent check of the support characterization in Theorem 3: cross-validation between the independently implemented vertex-indexed and block-indexed series computations, after the corresponding change of indexing. The eight *extended* rows additionally get the PARI/GP cross-check, and the numeric critical-pair consistency check against Theorem 10 or Proposition 11. This is not, and cannot be, a test of "all" finite Ω — there are infinitely many — but between the sizes, parities, and (for the mixed cases) the number of even lengths present, every qualitative combination the theorems distinguish is now exercised at least once: {5,6,7} (Table 1) is the |Ω|=3 example with exactly one even length, and {6,8,9} (extended) is the |Ω|=3 example with two. The extended eight are verification cases only; they are not additional theoretical results claimed in paper [2].

### What each cross-check verifies

| File | Independent route |
|---|---|
| `growth_rate_omega.py` | Estimates ρ_Ω from a plain ratio of consecutive series coefficients, instead of solving the criticality condition analytically. |
| `exhaustive_iso_omega.py` | For the smallest cases (1–2 blocks) of Ω={5,6}, enumerates the cactus graphs directly with `networkx` and counts isomorphism classes — no functional equation at all. |
| `exhaustive_iso_omega_extended.py` *(extended)* | The same direct-construction check, on Ω={7,11} instead. |
| `split_tree_omega.py` | Checks by brute force that the split-decomposition characterization (Theorem 1) still holds when two *different* cycle sizes from Ω={5,6} meet at the same cut vertex — the mixed-length case introduced by paper [2]. |
| `split_tree_omega_extended.py` *(extended)* | The same brute-force check, on Ω={7,11} instead. |
| `verify_dissymmetry_omega.py` | Recomputes the unrooted series from scratch using `sympy.Rational` and hand-written convolutions, instead of `fractions.Fraction`. Runs any Ω via `--omega`; used on Table 1's four Ω. |
| `verify_pari_omega.gp` | A third full implementation, in PARI/GP instead of Python, using PARI's native truncated power-series arithmetic, for Ω={5,6}. |
| `verify_extended_omega.py` *(extended)* | For eight Ω not in Table 1: Theorem 3's support check, agreement between the two independent Python series implementations, and consistency of the numerically located critical pair with Theorem 10 / Proposition 11. |
| `verify_pari_omega_extended.gp` *(extended)* | The same eight Ω through PARI/GP, cross-checked against `verify_extended_omega.py`'s output. |
| `verify_R4_recurrence.py` | Re-derives the rooted enumeration coefficients s_n from an independent recurrence obtained by logarithmic differentiation of Eq. (1) of the paper, built up one coefficient at a time (triangular, not circular), and checks exact agreement with the main fixed-point solver over a wide range of n, for Table 1's four Ω. |

### Independent cross-checks of the asymptotic law

Two further scripts check part (b) of Theorem 14 — the sub-exponential amplitude constant C_Ω of the asymptotic law — by a route independent of `critical_point_solver.py`'s own internal computation, for the four Ω of Table 1:

| File | Independent route |
|---|---|
| `amplitude_check_omega.py` | Computes C_Ω two independent ways — directly in the original coordinates, and via the reduced functional equation — and checks the two formulations agree, to a working precision of 50 significant digits. |
| `convergence_Rn_check.py` | Checks directly from the combinatorial recurrence for s_n (no closed form used) that R_n = s_n·ρ_Ω^n·n^(3/2) converges toward C_Ω, for all four Table 1 Ω, up to n=300. A numerical convergence check, not a proof of the convergence rate. |

This gives ten independent structural/enumeration cross-check scripts (the table above, including `verify_R4_recurrence.py`), plus these two additional asymptotic-amplitude checks, plus the two PARI/GP implementations.

## Tests (`tests/`)

**→ Run the whole suite with `python3 -m pytest tests/ -v`** (see Usage) — it locks in the results above against silent regressions. The suite currently comprises 14 tests across the eight files below.

- `test_regression_known_values.py` - solver output anchored against values already published in [1] (Theorem 2 and Table 3), covering both parities of the kernel.
- `test_cross_consistency.py` - agreement between the exact-series route and the direct numerical solver.
- `test_theorem_a_untabulated_omega.py` - paper [2]'s support-characterization theorem (Theorem 3), verified on Ω not appearing in its published Table 1.
- `test_untabulated_omega_extended.py` - the two fast checks of `verify_extended_omega.py` (Theorem 3, cross-implementation agreement), run automatically for eight Ω not appearing in Table 1.
- `test_random_omega_property.py` - Theorem 3's support characterization, on 30 finite sets Ω drawn at random (fixed seed, |Ω| ∈ {2,3,4}) rather than relying only on hand-selected examples.
- `test_boundary_single_length_omega.py` - the |Ω|=1 boundary case: checks that the general block-indexed series, applied to a singleton Ω, reproduces series independently published and verified on the OEIS (A398033, A397250, A397210) — not just the critical pair already anchored by `test_regression_known_values.py`, but the full series.
- `test_vertex_indexed_regression.py` - locks the vertex-indexed engine (`mgonal_cactus_series_omega.py`) against silent regressions.
- `test_oeis_data_regression.py` - locks the eight sequences prepared for OEIS submission (block-indexed) against silent regressions.

## Scope and limitations

The mathematical proofs are given in paper [2]; this repository does not replace them. As the paper itself states of its numerical tables and figure, these computations "are independent of the proofs of the stated theorems." The repository's role is to reproduce their numerical and combinatorial consequences, generate the tabulated sequences, and check the implementations against independent computational routes.

- The paper and this repository work in the regime *m* ≥ 5. The restriction comes from the split-decomposition framework used there (for *m* ≥ 5, the *m*-cycle is prime in the relevant sense; for *m* = 4 it is not); cycle lengths *m* = 3, 4 are outside its scope and are not treated here.
- The twelve-Ω coverage table above is a deliberately varied sample, not an exhaustive test: there are infinitely many finite Ω, and no finite set of cross-checks can cover them all. The sample was chosen to exercise every qualitative case the paper's theorems distinguish (size, parity, and — for mixed cases — the number of even lengths present), not to approximate exhaustiveness.
- The asymptotic-amplitude cross-checks (`amplitude_check_omega.py`, `convergence_Rn_check.py`) and the coefficient-recurrence cross-check (`verify_R4_recurrence.py`) currently cover only the four Ω of Table 1, not the eight extended Ω.
- `convergence_Rn_check.py` verifies convergence of R_n numerically, up to n=300; this is evidence for, not a proof of, the convergence rate, which the paper establishes analytically (Theorem 14(a)).
- The exhaustive direct-construction checks (`exhaustive_iso_omega.py`, `exhaustive_iso_omega_extended.py`, `split_tree_omega.py`, `split_tree_omega_extended.py`) are limited to the smallest block counts (1–2 blocks) by construction; beyond that, direct enumeration becomes intractable, which is exactly why the rest of the pipeline exists.

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
python3 python/verify_R4_recurrence.py
python3 python/amplitude_check_omega.py
python3 python/convergence_Rn_check.py
gp -q pari/verify_pari_omega.gp < /dev/null

# extended independent cross-checks (eight further, untabulated Omega)
python3 python/verify_extended_omega.py   # several minutes; not part of pytest
python3 python/exhaustive_iso_omega_extended.py
python3 python/split_tree_omega_extended.py
gp -q pari/verify_pari_omega_extended.gp < /dev/null
```

## Data availability

There are four mixed sets Ω covered here — {5,6}, {5,7}, {5,7,9}, {5,6,7}. For each one, Table 1 of paper [2] gives a summary of its growth parameters (ρ_Ω, τ_Ω, and whether an exact closed-form characterization of the critical value is available) — one row per Ω, four rows in total. Separately, for each of those same four Ω, the *full* rooted and unrooted term-by-term counts are computed and prepared for OEIS — two sequences per Ω, eight sequences in total. The four-row summary and the eight full sequences describe the same four cases at two different levels of detail, not eight different cases; the [companion page](docs/index.html) keeps them in two separate tables for that reason, rather than mixing a growth-parameter summary with raw term data in one table.

All eight sequences were computed by `mgonal_cactus_series_omega_blocks.py`, indexed by number of blocks per the convention this OEIS family already uses. These block-indexed rooted and unrooted sequences are complementary to the vertex-indexed series s_Ω(x) studied in the paper, not a restatement of it — see the *Method note* above for how the two indexings relate. 100 verified terms per sequence are provided in `bfiles/`, one per line in the standard OEIS b-file format (`index value` pairs, index starting at 1).

All eight have been reviewed and approved on OEIS, each carrying the data, a b-file, and a reference to this repository and to paper [2]:

| Ω | rooted | unrooted |
|---|---|---|
| {5,6} | [A399365](https://oeis.org/A399365) | [A399366](https://oeis.org/A399366) |
| {5,7} | [A399555](https://oeis.org/A399555) | [A399556](https://oeis.org/A399556) |
| {5,7,9} | [A397121](https://oeis.org/A397121) | [A399713](https://oeis.org/A399713) |
| {5,6,7} | [A399759](https://oeis.org/A399759) | [A399876](https://oeis.org/A399876) |

### Two exact results behind these sequences

- **A closed-form characterization in the all-odd case.** For the all-odd cases (Ω = {5,7} and Ω = {5,7,9}), the critical value τ_Ω is characterized by the exact closed-form equation of Theorem 10 of [2], rather than only by a numerical solution.
- **An exact description of where the sequence is nonzero.** For the vertex-indexed version of each sequence, the full set of nonzero positions is characterized (Theorem 3 of [2]) via the numerical semigroup generated by {m − 1 : m ∈ Ω} — reducing, for two generators, to the classical Frobenius coin problem.

## References

[1] Speyser, F. G. *Enumeration and Asymptotic Analysis of Strict Non-Plane m-Gonal Cactus Graphs via Split-Decomposition*, 2026. Preprint: DOI [10.5281/zenodo.21513753](https://zenodo.org/records/21513753).

[2] Speyser, F. G. & Vyatkina, K. *Enumeration and Asymptotic Analysis of Strict Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths*, 2026. Preprint: Zenodo. DOI: [10.5281/zenodo.22875493](https://doi.org/10.5281/zenodo.22875492).

## Citation

If you use this code, please cite it via its citable archive on Zenodo: DOI: [10.5281/zenodo.21854629](https://doi.org/10.5281/zenodo.21854629) (this DOI always resolves to the latest archived version of the code).

## Authors

Frédéric G. Speyser - Independent Researcher, Association Sciences & Coopération, France
ORCID: [0000-0002-1767-5325](https://orcid.org/0000-0002-1767-5325)

Kseniya Vyatkina - Independent Researcher, Association Sciences & Coopération, France
ORCID: [0009-0001-8849-6332](https://orcid.org/0009-0001-8849-6332)

## License

MIT (see `LICENSE`), for consistency with the other repositories.
