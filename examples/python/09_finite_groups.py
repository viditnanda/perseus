"""Homology and cohomology rings of finite groups, and persistence along a series of groups.

The input files are written from GAP by the exporters in tools/ (see examples/README.md); the three used
here are in examples/data. Run:  python3 examples/python/09_finite_groups.py
"""
import os

import perseus

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


def groups(hom):
    out = []
    for h in hom:
        parts = ([f"Z^{h['free']}" if h["free"] > 1 else "Z"] if h["free"] else []) + [f"Z/{d}" for d in h["torsion"]]
        out.append(f"H_{h['degree']} = " + (" + ".join(parts) or "0"))
    return ", ".join(out)


# 1. A solvable group from a polycyclic presentation: the generalised quaternion group Q_16
#    (gap> WritePerseusPcGroup(SmallGroup(16, 9), "q16_pcgroup.txt"); after Read("tools/export_pcgroup.g")).
r = perseus.pc_group(**perseus.read_pcgroup(os.path.join(DATA, "q16_pcgroup.txt")), maxdim=6, ring=2)
print(f"Q_16 (order {r['order']}):", groups(r["homology"]))
c = r["cohomology"]
print("  dim H^n(G; F_2), n = 0..6:", c["dims"])
print("  minimal ring generators (degree, index):", c["generators"])
print("  e.g. the square of the first degree-1 generator, in the basis of H^2:", c["products"][(1, 0, 1, 0)])
print("  pcgs used (relative orders):", r["pcgs_homology"]["relative_orders"], "-", r["pcgs_homology"]["choice"])

# 2. Any finite group, through its Sylow subgroups: the alternating group A_5
#    (gap> WritePerseusSylow(AlternatingGroup(5), "a5_sylow.txt"); after Read("tools/export_sylow.g")).
r = perseus.finite_group(**perseus.read_sylow(os.path.join(DATA, "a5_sylow.txt")), maxdim=6, ring=2)
print(f"\nA_5 (order {r['order']}, primes {r['primes']}):", groups(r["homology"]))
print("  dim H^n(A_5; F_2), n = 0..6:", r["cohomology"]["dims"], "generators by degree:", r["cohomology"]["generators_by_degree"])

# 3. Persistence along a series: the quotients of the dihedral group of order 64 by its lower central
#    series, D_64 -> D_32 -> D_16 -> D_8 -> C_2 x C_2, and the maps they induce on H_n(-; F_2)
#    (gap> WritePerseusQuotientSeries(DihedralGroup(64), LowerCentralSeries(DihedralGroup(64)), "...")).
r = perseus.pc_series(**perseus.read_pcseries(os.path.join(DATA, "d64_lower_central.txt")), maxdim=3, p=2)
print("\nD_64 -> D_32 -> D_16 -> D_8 -> C_2 x C_2 (quotients by the lower central series):")
for d in r["degrees"]:
    bars = ", ".join(f"[{i + 1}, {j + 1}] x {k}" for i, j, k in d["barcode"])
    print(f"  H_{d['degree']}(-; F_2) along the series: ranks of the composites", d["rank"][0], "| barcode", bars)
print("  (groups numbered 1..5; [i, j] x k: k classes born in G_i that survive through G_j)")
