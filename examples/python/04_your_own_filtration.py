"""Filtrations you build yourself: simplicial complexes and weighted graphs.

`simplices` takes (vertices, birth) pairs; faces that are not listed are added, born with their earliest
coface. `flag` takes a weighted graph and builds its clique complex. Run:
    python3 examples/python/04_your_own_filtration.py
"""
import os
from itertools import combinations

import numpy as np
import perseus

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")

# 1. The lower-star filtration of a height function on a triangulated torus: every simplex is born at the
#    largest height of its vertices. The barcode records the critical values: a minimum (H0), two saddles
#    (H1) and a maximum (H2), here at heights -3, -1, 1 and 3.
n = 8
idx = lambda i, j: (i % n) * n + (j % n)
height = {idx(i, j): (2 + np.cos(2 * np.pi * j / n)) * np.cos(2 * np.pi * i / n) for i in range(n) for j in range(n)}
tris = [t for i in range(n) for j in range(n)
        for t in ((idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)), (idx(i, j), idx(i, j + 1), idx(i + 1, j + 1)))]
simplices = {tuple(sorted(f)) for t in tris for k in (1, 2, 3) for f in combinations(t, k)}
S = [(s, max(height[v] for v in s)) for s in simplices]
r = perseus.simplices(S)
print("torus, lower-star filtration of the height:")
for k, B in enumerate(r["barcode"]):
    print(f"  H{k}:", np.round(B, 3).tolist())

# 2. The clique (flag) filtration of a weighted graph: a 4-cycle a-b-c-d whose chord a-c appears later,
#    with vertex d born late. Vertices are integers here; examples/data/square_flag.txt has the same graph
#    with labels for the command line.
a, b, c, d = range(4)
edges = [(a, b, 1.0), (b, c, 1.0), (c, d, 1.0), (d, a, 1.0), (a, c, 2.0)]
r = perseus.flag(edges, vertex_births=[0, 0, 0, 0.5], maxdim=2)
print("\nweighted 4-cycle:", {k: np.round(B, 3).tolist() for k, B in enumerate(r["barcode"])})
print("  the loop is born at 1 and filled at 2, when the chord makes two triangles")

# 3. The boundary of a tetrahedron (a 2-sphere) whose faces arrive one by one. Unlisted edges and vertices
#    come with their first triangle, so the complex is a growing disc until the last face closes the sphere.
faces = list(combinations(range(4), 3))
r = perseus.simplices([(f, float(t)) for t, f in enumerate(faces)])
print("\nboundary of a tetrahedron, faces at times 0..3:", {k: B.tolist() for k, B in enumerate(r["barcode"])})
