"""Magnitude homology of graphs over Z, torsion included (Hepworth-Willerton).

MH_{k,l}(G; Z) is computed exactly for lengths l <= lmax; `pairs` restricts to chosen endpoints (s, t).
Run:  python3 examples/python/07_magnitude.py
"""
import os

import perseus

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


def show(table):
    for (k, l), g in sorted(table.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        parts = ([f"Z^{g['free']}"] if g["free"] else []) + [f"(Z/{d})^{m}" if m > 1 else f"Z/{d}" for d, m in g["torsion"].items()]
        if parts:
            print(f"  MH_{{{k},{l}}} = " + " + ".join(parts))


# 1. The 7-cycle: free groups only.
C7 = [(i, (i + 1) % 7) for i in range(7)]
print("cycle graph C_7, l <= 5:")
show(perseus.magnitude(C7, lmax=5)["table"])

# 2. Torsion. The Kaneta-Yoshinaga graph of a regular CW complex M is the Hasse diagram of its face poset
#    with a bottom and a top added; the reduced homology of M embeds in its magnitude homology between
#    bottom and top. For the 6-vertex RP^2, the Z/2 of H_1(RP^2) appears as MH_{3,4}.
E = [tuple(map(int, line.split())) for line in open(os.path.join(DATA, "ky_rp2.txt")) if line.strip() and not line.startswith("#")]
n = 1 + max(max(e) for e in E)
bottom, top = n - 2, n - 1
print(f"\nKaneta-Yoshinaga graph of RP^2 ({n} vertices), between bottom and top, l <= 6:")
show(perseus.magnitude(E, lmax=6, pairs=[(bottom, top)])["table"])

# Directed graphs: directed=True.
print("\ndirected 4-cycle, l <= 4:")
show(perseus.magnitude([(0, 1), (1, 2), (2, 3), (3, 0)], directed=True, lmax=4)["table"])
