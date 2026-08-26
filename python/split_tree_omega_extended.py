"""
split_tree_omega_extended.py

The same brute-force split-decomposition check as split_tree_omega.py,
run on a heptagon and a hendecagon (Omega = {7, 11}) sharing a cut
vertex, instead of a pentagon and a hexagon (Omega = {5, 6}). This
confirms the mixed-size positive/negative test pattern is not somehow
specific to the {5,6} sizes checked in split_tree_omega.py.

Reference: F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths" (in
preparation for submission to the Journal of Integer Sequences, 2026).

Auteur : Frederic G. Speyser
Run: python3 split_tree_omega_extended.py
"""
from itertools import combinations


def neighbors(edges, v):
    return {b for a, b in edges if a == v} | {a for a, b in edges if b == v}


def find_split(vertices, edges):
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


def main():
    print("=" * 70)
    print("POSITIVE TEST: a heptagon and a hendecagon sharing one cut vertex")
    print("(Omega = {7, 11}, untabulated)")
    print("=" * 70)
    heptagon = list(range(7))
    hendecagon = [0] + list(range(20, 30))   # 11 vertices total, shares vertex 0
    edges = cycle_edges(heptagon) + cycle_edges(hendecagon)
    vertices = set(heptagon) | set(hendecagon)
    split = find_split(vertices, edges)
    print(f"Split found: {split is not None}")
    ok_hept, msg_hept = is_clean_m_cycle(
        heptagon, block_edges(vertices, edges, heptagon), 7)
    ok_hend, msg_hend = is_clean_m_cycle(
        hendecagon, block_edges(vertices, edges, hendecagon), 11)
    print(f"  Heptagon block (by construction) is a clean C_7: "
          f"{ok_hept} ({msg_hept})")
    print(f"  Hendecagon block (by construction) is a clean C_11: "
          f"{ok_hend} ({msg_hend})")
    positive_ok = ok_hept and ok_hend
    print(f"  ==> Theorem 1, condition (a) holds for a MIXED-size pair: "
          f"{positive_ok}")

    print()
    print("=" * 70)
    print("NEGATIVE TEST #1: same graph + a CHORD added inside the heptagon")
    print("=" * 70)
    edges_chord = edges + [(1, 4)]
    side_with_chord = block_edges(vertices, edges_chord, heptagon)
    ok_chord, msg_chord = is_clean_m_cycle(heptagon, side_with_chord, 7)
    print(f"  Heptagon with a chord is a clean C_7: {ok_chord} ({msg_chord})")
    print(f"  ==> Theorem 1, condition (a) VIOLATED, as expected: {not ok_chord}")

    print()
    print("=" * 70)
    print("NEGATIVE TEST #2: same graph + a BRIDGE between the two differently-")
    print("sized blocks")
    print("=" * 70)
    edges_bridge = edges + [(3, 25)]
    split_bridge = find_split(vertices, edges_bridge)
    print(f"  Split found despite the bridge: {split_bridge is not None}")
    print("  ==> the bridge destroys the 'strict' structure regardless of the")
    print("      two block sizes being different: this graph is no longer in")
    print("      the class Theorem 1 covers to begin with -- expected.")

    print()
    print("=" * 70)
    print("NEGATIVE TEST #3: a cycle length NOT in Omega={7,11} attached at the")
    print("shared vertex (a hexagon instead of a hendecagon)")
    print("=" * 70)
    hexagon = [0] + list(range(40, 45))
    edges_wrong = cycle_edges(heptagon) + cycle_edges(hexagon)
    vertices_wrong = set(heptagon) | set(hexagon)
    ok_w7, msg_w7 = is_clean_m_cycle(
        hexagon, block_edges(vertices_wrong, edges_wrong, hexagon), 7)
    ok_w11, msg_w11 = is_clean_m_cycle(
        hexagon, block_edges(vertices_wrong, edges_wrong, hexagon), 11)
    print(f"  Hexagon block is a clean C_7: {ok_w7} ({msg_w7})")
    print(f"  Hexagon block is a clean C_11: {ok_w11} ({msg_w11})")
    print(f"  ==> Neither matches -- correctly rejected as outside Omega={{7,11}}: "
          f"{not ok_w7 and not ok_w11}")

    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    all_ok = (positive_ok and not ok_chord and
              split_bridge is None and not ok_w7 and not ok_w11)
    print(f"All checks behaved as expected: {all_ok}")
    return all_ok


if __name__ == "__main__":
    ok = main()
    import sys
    sys.exit(0 if ok else 1)
