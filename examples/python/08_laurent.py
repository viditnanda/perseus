"""Chain complexes over the Laurent ring Z[t, t^-1]: infinite cyclic covers, Alexander polynomials,
Novikov homology.

Each cell is (dim, birth, [(face, polynomial), ...]) with polynomials in t as strings.
Run:  python3 examples/python/08_laurent.py
"""
import perseus

# The trefoil complement from the presentation <x, y | xyx = yxy>, lifted to its infinite cyclic cover by
# Fox calculus: one 0-cell, two 1-cells (x and y, each with boundary (t - 1) v), one 2-cell r with
# boundary (1 - t + t^2) x - (1 - t + t^2) y. Its H_1 over Q[t, t^-1] is Q[t, t^-1]/(Alexander polynomial).
trefoil = [
    (0, 0, []),
    (1, 0, [(0, "t - 1")]),
    (1, 0, [(0, "t - 1")]),
    (2, 0, [(1, "t^2 - t + 1"), (2, "-t^2 + t - 1")]),
]
r = perseus.laurent(trefoil)
print("trefoil, homology over", r["homology_over"])
for h in r["homology"]:
    torsion = " + ".join(f"L/({f})" for f in h["torsion"])
    print(f"  H_{h['degree']} =", " + ".join(p for p in ([f"L^{h['free']}"] if h["free"] else []) + ([torsion] if torsion else [])) or "0")
print("  (L = Q[t, t^-1]; H_1 = L/(t^2 - t + 1): the Alexander polynomial of the trefoil)")

# Over the field Q(t) (equivalently the Novikov field Q((t))) everything is torsion, so the barcode is empty.
# The irreducibles f for which the residue field Q[t]/(f) gives a different barcode are exactly the
# factors of the torsion: here t - 1 and the Alexander polynomial.
print("barcode over", r["barcode_over"], {k: B.tolist() for k, B in enumerate(r["barcode"])})
for b in r["bad_irreducibles"]:
    print(f"  bad irreducible {b['polynomial']}:", {k: B.tolist() for k, B in enumerate(b["barcode"])})
print("bad primes for F_p(t):", sorted(r["bad_primes"]) or "none")

# Integral content: Z[t, t^-1]/(2, t + 1) as H_0 of a vertex v with two edges, of boundary 2 v and (t + 1) v.
# Over Q[t, t^-1] the 2 is a unit and H_0 = 0; over F_2[t, t^-1] it is not, and H_0 = F_2[t, t^-1]/(t + 1).
cells = [(0, 0, []), (1, 0, [(0, "2")]), (1, 0, [(0, "t + 1")])]
r = perseus.laurent(cells)
print("\nZ[t, t^-1]/(2, t + 1): H_0 over Q[t, t^-1] has free rank", r["homology"][0]["free"],
      "and torsion", r["homology"][0]["torsion"] or "none")
for p, hom in r["homology_bad_primes"].items():
    print(f"  but over F_{p}[t, t^-1]:", ", ".join(
        f"H_{h['degree']} = " + (" + ".join(([f"L^{h['free']}"] if h["free"] else []) + [f"L/({f})" for f in h["torsion"]]) or "0")
        for h in hom), f"(L = F_{p}[t, t^-1])")
