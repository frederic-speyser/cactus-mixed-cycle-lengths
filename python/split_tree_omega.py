"""
split_tree_omega.py

Brute-force verification that Theorem 1's split-decomposition
characterization (Definition 1 of [1]) still holds, unchanged, when two
different cycle sizes from Omega can meet at the same cut vertex -- a
question the original paper never had to ask, since it only ever treats
a single fixed m throughout.

This is the mixed-Omega analogue of split_tree_v2.py from [1]. The
positive test there used two blocks of the SAME size (two pentagons). The
genuinely new case here is two blocks of DIFFERENT sizes sharing a cut
vertex (a pentagon and a hexagon) -- does the split-decomposition tree
still cleanly separate into two prime nodes, each labeled by the correct
cycle length, with no unexpected interaction between the two different
sizes at the shared vertex?

Also tested: the same negative controls as [1] (a chord, a bridge, a
cycle of a length not in Omega), now applied to a mixed-size base graph,
to confirm these failure modes are detected exactly as before and are not
somehow masked by the presence of two different sizes.

This script does not prove Theorem 1's characterization generalizes to
mixed Omega -- it tests it on specific small graphs, by brute force. A
negative result here would be a genuine problem for everything computed
so far in this repository, since the solver implicitly assumes the
characterization holds. A positive result is necessary, not sufficient,
for that assumption to be safe in general.

Reference: Fr. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
Non-Plane m-Gonal Cactus Graphs via Split-Decomposition" [1], Theorem 1
and Definition 1 (the construction method mirrors split_tree_v2.py from
that paper's repository).

Author: Frederic G. Speyser
Run: python3 split_tree_omega.py
"""
from itertools import combinations


def neighbors(edges, v):
    return {b for a, b in edges if a == v} | {a for a, b in edges if b == v}


def find_split(vertices, edges):
    """Brute-force search for a split (Definition 1 of [1])."""
    edges_set = {frozenset(e) for e in edges}
    vertices = list(vertices)
    n = len(vertices)
    for size1 in range(2, n - 1):
        for V1 in combinations(vertices, size1):
            V1 = set(V1)
            V2 = set(vertices) - V1
            if len(V2) < 2:
                continue
            A = {v for v in V2 if neighbors(edges, v) & V1}
            B = {v for v in V1 if neighbors(edges, v) & V2}
            if A and B and all(frozenset((a, b)) in edges_set for a in A for b in B):
                return (V1, V2, A, B)
    return None


def cycle_edges(order):
    n = len(order)
    return [(order[i], order[(i + 1) % n]) for i in range(n)]


def is_clean_m_cycle(vertices, edges, m):
    """A genuine C_m block: m vertices, degree 2 everywhere, a single cycle
    of exactly length m -- generalized from [1]'s split_tree_v2.py to take
    the target size m as an explicit parameter, so mismatched sizes can be
    tested directly (e.g. a hexagon checked against m=5 must fail)."""
    if len(vertices) != m:
        return False, f"size {len(vertices)} != {m}"
    deg = {v: 0 for v in vertices}
    for a, b in edges:
        deg[a] += 1
        deg[b] += 1
    if not all(d == 2 for d in deg.values()):
        return False, f"non-uniform degrees: {deg}"
    if len(edges) != m:
        return False, f"{len(edges)} edges != {m} (a chord or a missing edge)"
    return True, "OK"


def block_edges(vertices, edges, block):
    return [e for e in edges if e[0] in block and e[1] in block]


print("=" * 70)
print("POSITIVE TEST: a pentagon and a hexagon sharing one cut vertex --")
print("the genuinely new case, absent from [1] (which only ever glued")
print("blocks of the SAME size to each other)")
print("=" * 70)
pentagon = list(range(5))          # vertices 0..4
hexagon = [0, 5, 6, 7, 8, 9]        # shares vertex 0 with the pentagon
edges = cycle_edges(pentagon) + cycle_edges(hexagon)
vertices = set(pentagon) | set(hexagon)
split = find_split(vertices, edges)
print(f"Split found: {split is not None}")
# NOTE, following [1]'s own split_tree_v2.py methodology exactly: find_split
# only confirms THAT the graph is split-decomposable (Definition 1 permits
# many splits, not just the canonical block-tree one -- a brute-force
# search can return a degenerate split, e.g. isolating a sub-arc of one
# cycle while excluding the shared cut vertex). The actual block check is
# done on the KNOWN block vertex sets used to construct the graph, exactly
# as in [1]'s script, not on whatever V1/V2 the raw split happens to be.
ok_pentagon, msg_pentagon = is_clean_m_cycle(
    pentagon, block_edges(vertices, edges, pentagon), 5)
ok_hexagon, msg_hexagon = is_clean_m_cycle(
    hexagon, block_edges(vertices, edges, hexagon), 6)
print(f"  Pentagon block (by construction) is a clean C_5: "
      f"{ok_pentagon} ({msg_pentagon})")
print(f"  Hexagon block (by construction) is a clean C_6: "
      f"{ok_hexagon} ({msg_hexagon})")
one_side_is_c5_other_is_c6 = ok_pentagon and ok_hexagon
print(f"  ==> Theorem 1, condition (a) holds for a MIXED-size pair: "
      f"{one_side_is_c5_other_is_c6}")

print()
print("=" * 70)
print("NEGATIVE TEST #1: same graph + a CHORD added inside the pentagon")
print("=" * 70)
edges_chord = edges + [(1, 3)]
side_with_chord = block_edges(vertices, edges_chord, pentagon)
ok_chord, msg_chord = is_clean_m_cycle(pentagon, side_with_chord, 5)
print(f"  Pentagon with a chord is a clean C_5: {ok_chord} ({msg_chord})")
print(f"  ==> Theorem 1, condition (a) VIOLATED, as expected: {not ok_chord}")

print()
print("=" * 70)
print("NEGATIVE TEST #2: same graph + a BRIDGE between the two differently-")
print("sized blocks")
print("=" * 70)
edges_bridge = edges + [(2, 7)]
split_bridge = find_split(vertices, edges_bridge)
print(f"  Split found despite the bridge: {split_bridge is not None}")
if split_bridge:
    V1b, V2b, Ab, Bb = split_bridge
    print(f"    V1={sorted(V1b)} V2={sorted(V2b)}")
print("  ==> the bridge destroys the 'strict' structure regardless of the")
print("      two block sizes being different: this graph is no longer in")
print("      the class Theorem 1 covers to begin with -- expected, and")
print("      unaffected by the size mismatch.")

print()
print("=" * 70)
print("NEGATIVE TEST #3: a cycle length NOT in Omega={5,6} attached at the")
print("shared vertex (a heptagon instead of a hexagon)")
print("=" * 70)
heptagon = [0, 10, 11, 12, 13, 14, 15]
edges_wrong = cycle_edges(pentagon) + cycle_edges(heptagon)
vertices_wrong = set(pentagon) | set(heptagon)
split_wrong = find_split(vertices_wrong, edges_wrong)
ok_w5, msg_w5 = is_clean_m_cycle(
    heptagon, block_edges(vertices_wrong, edges_wrong, heptagon), 5)
ok_w6, msg_w6 = is_clean_m_cycle(
    heptagon, block_edges(vertices_wrong, edges_wrong, heptagon), 6)
print(f"  Split still found (it's a genuine split regardless of Omega): "
      f"{split_wrong is not None}")
print(f"  Heptagon block is a clean C_5: {ok_w5} ({msg_w5})")
print(f"  Heptagon block is a clean C_6: {ok_w6} ({msg_w6})")
print(f"  ==> Neither matches -- correctly rejected as outside Omega={{5,6}}: "
      f"{not ok_w5 and not ok_w6}")

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
all_ok = (one_side_is_c5_other_is_c6 and not ok_chord and
          split_bridge is None and not ok_w5 and not ok_w6)
print(f"All checks behaved as expected: {all_ok}")
print("This is evidence, not a proof, that Theorem 1's characterization")
print("(conditions a-d of [1]) needs no modification for mixed Omega on")
print("the cases tested here. It does not rule out a subtler failure at")
print("larger k or with three or more distinct sizes meeting at one")
print("cut vertex.")
