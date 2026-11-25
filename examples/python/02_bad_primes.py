"""Coefficients matter: the barcode over Q, the bad primes, and what changes over F_p.

A prime p is *bad* when the barcode over F_p differs from the one over Q. Perseus computes the barcode
over Q, finds every bad prime, proves the list is complete, and says which intervals change.
Run:  python3 examples/python/02_bad_primes.py
"""
import os

import numpy as np
import perseus

show = lambda a: np.round(np.asarray(a, dtype=float), 3).tolist()   # intervals, rounded for printing
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")

# 1. The Moore space M(Z/12, 1): a vertex, a loop at time 1, and a disc attached along the loop 12 times
#    at time 2. Each cell is (dimension, birth, [(face, coefficient), ...]) with faces indexed by position.
moore = [(0, 0.0, []), (1, 1.0, []), (2, 2.0, [(1, 12)])]
r = perseus.cells(moore, fp=[2, 3, 5])
print("M(Z/12, 1)")
print("  over Q, H1:", r["barcode"][1].tolist(), "  (the disc kills the loop)")
print("  bad primes:", perseus.bad_primes(r), "  (the primes dividing 12)")
for p in (2, 3, 5):
    print(f"  over F_{p}, H1:", r["fields"][p][1].tolist())

# 2. RP^2 from its 6-vertex triangulation. Over Q it looks like a point; over F_2 it has H1 and H2.
S = []
for line in open(os.path.join(DATA, "rp2_simplices.txt")):
    if line.strip() and not line.startswith("#"):
        *v, birth = line.split()
        S.append((tuple(int(x) for x in v), float(birth)))
r = perseus.simplices(S)
print("\nRP^2 (6 vertices): bad primes", perseus.bad_primes(r))
print("  intervals only over F_2, by degree:", [show(a) for a in r["certificate"]["bad_primes"][2]["only_over_Fp"]])

# 3. RP^2 sampled as a point cloud (the Veronese surface in R^6). The longest H1 interval over F_2 is the
#    nonorientable loop; over Q the same class dies early, killed by twice a 2-cell.
X = np.loadtxt(os.path.join(DATA, "rp2_points.txt"))
r = perseus.rips(points=X, maxdim=2)
bp = r["certificate"]["bad_primes"]
print(f"\nRP^2, {len(X)} sample points: bad primes {perseus.bad_primes(r)}; proved: {r['certificate']['proof']}")
print("  the H1 interval that differs, over Q :", show(bp[2]["only_over_Q"][1]))
print("  the same H1 interval over F_2        :", show(bp[2]["only_over_Fp"][1]))
print("  H2 intervals that exist only over F_2:", show(bp[2]["only_over_Fp"][2]))

# The non-unit pivots of the rational reduction are where the torsion is detected.
for piv in r["certificate"]["nonunit_pivots"]:
    print(f"  pivot {piv['value']} in degree {piv['degree']} at [{piv['birth']:.3f}, {piv['death']:.3f})")
