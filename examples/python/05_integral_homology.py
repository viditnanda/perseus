"""Integral information: homology over Z, images of maps, and the torsion along a filtration.

Over Z there is no barcode, but Perseus reports what can be said exactly: H_*(K; Z) with its torsion, the
image of H_*(K_s; Z) -> H_*(K_t; Z), and, for every value t, the torsion subgroup of H_*(K_t; Z).
Run:  python3 examples/python/05_integral_homology.py
"""
import os

import perseus

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


def groups(hom):
    """Format a list of {"free": r, "torsion": [d1, d2, ...]} by degree as Z^r + Z/d1 + ..."""
    out = []
    for k, h in enumerate(hom):
        parts = ([f"Z^{h['free']}" if h["free"] > 1 else "Z"] if h["free"] else []) + [f"Z/{d}" for d in h["torsion"]]
        out.append(f"H_{k} = " + (" + ".join(parts) or "0"))
    return ", ".join(out)


# 1. A two-stage filtration: RP^1 (a vertex and a loop) at time 0, then RP^2 (a disc attached by degree 2)
#    at time 1. The generator of H_1(RP^1) = Z maps onto H_1(RP^2) = Z/2.
cells = [(0, 0.0, []), (1, 0.0, []), (2, 1.0, [(1, 2)])]
r = perseus.cells(cells, integral=True, image=(0, 1))
I = r["integral"]
print("RP^1 -> RP^2")
print("  homology of the whole complex:", groups(I["homology"]))
print("  image of H_*(K_0; Z) -> H_*(K_1; Z):", groups(I["image"]["groups"]))

# 2. A random 2-complex (Linial-Meshulam process on 40 vertices): all edges at time 1, then triangles one
#    at a time. As H_1 dies it passes through torsion, which --torsion tracks exactly at every step.
S = []
for line in open(os.path.join(DATA, "lm40.txt")):
    if line.strip() and not line.startswith("#"):
        *v, t = line.split()
        S.append((tuple(int(x) for x in v), float(t)))
r = perseus.simplices(S, torsion=True)
print("\nLinial-Meshulam, 40 vertices, 700 triangles: bad primes", perseus.bad_primes(r))
for step in r["torsion"]["degrees"][1]:
    factors = step["invariant_factors"]
    print(f"  t = {step['value']:g}: tors H_1(K_t; Z) =", " + ".join(f"Z/{d}" for d in factors) or "0")
print("  the bad primes 3, 5, 13 are exactly the primes of this torsion (975 = 3 * 5^2 * 13)")
