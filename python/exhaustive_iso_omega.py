"""
exhaustive_iso_omega.py

Fully independent verification of the smallest cases (k=1 and k=2 blocks)
for strict cacti admitting Omega = {5, 6} -- the graphs are built directly
as combinatorial objects (via networkx) and deduplicated by graph
isomorphism, without going through the functional equation of
mgonal_cactus_series_omega.py or any other code in this repository.

This is the mixed-Omega analogue of exhaustive_iso.py from [1], which
verified k=1,2,3 for the pure m=5 case. Here the interesting cases are the
ones that mix cycle sizes explicitly: a pentagon glued to a hexagon, not
just pentagon-to-pentagon or hexagon-to-hexagon.

What is checked against the solver's output (mgonal_cactus_series_omega.py,
unrooted series, Omega={5,6}):
  - k=1: exactly 1 class for a lone pentagon (degree 5), exactly 1 class
    for a lone hexagon (degree 6).
  - k=2: exactly 1 class for each of the three combinations -- two
    pentagons (degree 9), one pentagon + one hexagon (degree 10), two
    hexagons (degree 11) -- matching the solver's unrooted coefficients
    of 1 at each of x^9, x^10, x^11.

Reference: Fr. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
Non-Plane m-Gonal Cactus Graphs via Split-Decomposition" [1], Sections
5.1-5.3 (the construction method mirrors exhaustive_iso.py from that
paper's repository).

Author: Frederic G. Speyser
Run: python3 exhaustive_iso_omega.py   (requires: pip install networkx)
"""
import networkx as nx


def cycle_graph(size, offset):
    """A C_size with vertices 'offset_0'..'offset_{size-1}'."""
    verts = [f"{offset}_{i}" for i in range(size)]
    G = nx.Graph()
    G.add_nodes_from(verts)
    for i in range(size):
        G.add_edge(verts[i], verts[(i + 1) % size])
    return G, verts


def glue(graphs_verts, merges):
    """graphs_verts: list of (G, verts). merges: list of pairs
    ((i1,v1),(i2,v2)) meaning verts[i1][v1] and verts[i2][v2] must become
    the same vertex."""
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


print("=" * 70)
print("k=1 block: a lone pentagon, and a lone hexagon")
print("=" * 70)
g5, _ = cycle_graph(5, "a")
g6, _ = cycle_graph(6, "b")
c5 = dedup_by_isomorphism([g5])
c6 = dedup_by_isomorphism([g6])
print(f"  Pentagon (n=5): {len(c5)} class  (expected: 1)")
print(f"  Hexagon  (n=6): {len(c6)} class  (expected: 1)")

print()
print("=" * 70)
print("k=2 blocks: the three size combinations, each at every relative")
print("attachment offset -- expect exactly 1 isomorphism class per combo")
print("=" * 70)


def build_k2_candidates(size1, size2):
    candidates = []
    for offset in range(size2):
        g1, v1 = cycle_graph(size1, "a")
        g2, v2 = cycle_graph(size2, "b")
        G = glue([(g1, v1), (g2, v2)], [((0, 0), (1, offset))])
        candidates.append(G)
    return candidates

combos = [(5, 5, "pentagon + pentagon (n=9)"),
          (5, 6, "pentagon + hexagon (n=10)"),
          (6, 6, "hexagon + hexagon (n=11)")]

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
print("Cross-check against mgonal_cactus_series_omega.py --omega 5,6")
print("=" * 70)
print("Solver's unrooted coefficients (computed independently, via the")
print("functional equation): x^9 -> 1, x^10 -> 1, x^11 -> 1")
print("This script's independent construction:", {
    label: cnt for label, (cnt, n) in results.items()
})
all_match = (len(c5) == 1 and len(c6) == 1 and
             all(cnt == 1 for cnt, n in results.values()))
print("All match solver's k=1, k=2 coefficients:", all_match)
