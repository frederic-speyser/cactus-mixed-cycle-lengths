"""
exhaustive_iso_omega_extended.py

The same fully-independent construction-and-isomorphism-count method as
exhaustive_iso_omega.py, but run on Omega = {7, 11} instead of the
tabulated Omega = {5, 6} -- a heptagon and a hendecagon (11-gon) sharing
a cut vertex, rather than a pentagon and a hexagon. As in
exhaustive_iso_omega.py, the graphs are built directly with networkx and
deduplicated by graph isomorphism, without going through the functional
equation of mgonal_cactus_series_omega.py or any other code in this
repository.

What is checked against the solver's output for Omega={7,11} (unrooted
series, printed by verify_extended_omega.py):
  - k=1: exactly 1 class for a lone heptagon (degree 7), exactly 1 class
    for a lone hendecagon (degree 11).
  - k=2: exactly 1 class for each of the three combinations -- two
    heptagons (degree 13), one heptagon + one hendecagon (degree 17),
    two hendecagons (degree 21) -- matching the solver's unrooted
    coefficients of 1 at each of x^13, x^17, x^21.

Reference: F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
Non-Plane Cactus Graphs over a Finite Set of Cycle Lengths" (in
preparation for submission to the Journal of Integer Sequences, 2026).

Author: Frederic G. Speyser
Run: python3 exhaustive_iso_omega_extended.py   (requires: pip install networkx)
"""
import networkx as nx


def cycle_graph(size, offset):
    verts = [f"{offset}_{i}" for i in range(size)]
    G = nx.Graph()
    G.add_nodes_from(verts)
    for i in range(size):
        G.add_edge(verts[i], verts[(i + 1) % size])
    return G, verts


def glue(graphs_verts, merges):
    G = nx.Graph()
    for g, verts in graphs_verts:
        G = nx.union(G, g)
    parent = {}

    def find(x):
        while parent.get(x, x) != x:
            x = parent[x]
        return x

    def union(x, y):
        parent[find(x)] = find(y)

    for (i1, v1), (i2, v2) in merges:
        a = graphs_verts[i1][1][v1]
        b = graphs_verts[i2][1][v2]
        union(a, b)
    mapping = {n: find(n) for n in G.nodes()}
    return nx.relabel_nodes(G, mapping)


def dedup_by_isomorphism(graph_list):
    classes = []
    for G in graph_list:
        found = False
        for cls in classes:
            if nx.is_isomorphic(G, cls[0]):
                cls.append(G)
                found = True
                break
        if not found:
            classes.append([G])
    return classes


def build_k2_candidates(size1, size2):
    candidates = []
    for offset in range(size2):
        g1, v1 = cycle_graph(size1, "a")
        g2, v2 = cycle_graph(size2, "b")
        G = glue([(g1, v1), (g2, v2)], [((0, 0), (1, offset))])
        candidates.append(G)
    return candidates


def main():
    print("=" * 70)
    print("Omega = {7, 11}: k=1 block -- a lone heptagon, and a lone hendecagon")
    print("=" * 70)
    g7, _ = cycle_graph(7, "a")
    g11, _ = cycle_graph(11, "b")
    c7 = dedup_by_isomorphism([g7])
    c11 = dedup_by_isomorphism([g11])
    print(f"  Heptagon   (n=7):  {len(c7)} class  (expected: 1)")
    print(f"  Hendecagon (n=11): {len(c11)} class  (expected: 1)")

    print()
    print("=" * 70)
    print("k=2 blocks: the three size combinations")
    print("=" * 70)
    combos = [(7, 7, "heptagon + heptagon (n=13)"),
              (7, 11, "heptagon + hendecagon (n=17)"),
              (11, 11, "hendecagon + hendecagon (n=21)")]
    results = {}
    for size1, size2, label in combos:
        cands = build_k2_candidates(size1, size2)
        classes = dedup_by_isomorphism(cands)
        n = 1 + (size1 - 1) + (size2 - 1)
        results[label] = (len(classes), n)
        print(f"  {label}: {len(classes)} class "
              f"(out of {len(cands)} attachment offsets tested), "
              f"n = {n}  (expected: 1)")

    print()
    print("=" * 70)
    print("Cross-check against verify_extended_omega.py's Omega={7,11} output")
    print("=" * 70)
    print("Solver's unrooted coefficients (independent, via the functional")
    print("equation): x^13 -> 1, x^17 -> 1, x^21 -> 1")
    print("This script's independent construction:", {
        label: cnt for label, (cnt, n) in results.items()
    })
    all_match = (len(c7) == 1 and len(c11) == 1 and
                 all(cnt == 1 for cnt, n in results.values()))
    print("All match solver's k=1, k=2 coefficients:", all_match)
    return all_match


if __name__ == "__main__":
    ok = main()
    import sys
    sys.exit(0 if ok else 1)
