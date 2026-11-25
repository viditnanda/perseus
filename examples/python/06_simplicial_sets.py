"""Simplicial sets, nerves of categories, and classifying spaces.

A simplicial set may glue a simplex to itself, so small models suffice: RP^2 needs one simplex in each
dimension. Run:  python3 examples/python/06_simplicial_sets.py
"""
from itertools import permutations

import perseus


def groups(hom):
    out = []
    for k, h in enumerate(hom):
        parts = ([f"Z^{h['free']}" if h["free"] > 1 else "Z"] if h["free"] else []) + [f"Z/{d}" for d in h["torsion"]]
        out.append(f"H_{k} = " + (" + ".join(parts) or "0"))
    return ", ".join(out)


# 1. RP^2 as a simplicial set: a vertex v (time 0), an edge e from v to v (time 0), and a triangle f with
#    faces d_0 f = e, d_1 f = s_0 v (the degenerate edge at v), d_2 f = e (time 1). Each simplex is
#    (dim, birth, [d_0, ..., d_dim]); a face is an earlier simplex y, or (y, [k, ...]) for s_k ... y.
rp2 = [(0, 0, []), (1, 0, [0, 0]), (2, 1, [1, (0, [0]), 1])]
r = perseus.simplicial_set(rp2, integral=True, image=(0, 1))
print("RP^2 as a simplicial set:", groups(r["integral"]["homology"]))
print("  bad primes:", perseus.bad_primes(r), "| image of H_*(circle) -> H_*(RP^2):", groups(r["integral"]["image"]["groups"]))

# 2. The nerve of a category: the face poset of a triangle's boundary (three vertices below three edges).
#    Its nerve is the barycentric subdivision, a circle. With weights, it is a filtration: the inclusions
#    into edge i are born at time i, so the loop closes at time 2.
objects = ["a", "b", "c", "ab", "bc", "ca"]
src, tgt, w = [], [], []
for o in range(len(objects)):                     # identities first: morphism o is the identity of object o
    src.append(o); tgt.append(o); w.append(0.0)
for t, e in enumerate(["ab", "bc", "ca"]):         # the inclusions vertex < edge
    for v in e:
        src.append(objects.index(v)); tgt.append(objects.index(e)); w.append(float(t))
comp = {}                                         # g o f for composable non-identities: there are none in a poset of height 1
r = perseus.nerve(src, tgt, list(range(len(objects))), comp, weights=w, maxdim=2, integral=True)
print("\nnerve of the face poset of a triangle's boundary:", groups(r["integral"]["homology"]))
print("  barcode:", {k: B.tolist() for k, B in enumerate(r["barcode"])})

# 3. Classifying spaces: the bar construction BG of a finite group (here the cyclic group Z/3 and the
#    symmetric group S_3) gives the group homology H_*(G; Z). It has (|G| - 1)^n simplices in degree n;
#    for solvable groups, pc_group (example 09) is far cheaper.
Z3 = [[(a + b) % 3 for b in range(3)] for a in range(3)]
print("\nH_*(Z/3; Z) from BZ/3:", groups(perseus.bar_construction(Z3, identity=0, maxdim=4)["integral"]["homology"]))
S3 = list(permutations(range(3)))
mult = [[S3.index(tuple(p[q[i]] for i in range(3))) for q in S3] for p in S3]
print("H_*(S_3; Z) from BS_3:", groups(perseus.bar_construction(mult, identity=S3.index((0, 1, 2)), maxdim=4)["integral"]["homology"]))
