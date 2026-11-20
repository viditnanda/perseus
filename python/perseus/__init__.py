"""Perseus 6: persistent homology over Q with certified bad primes, or over any prime field.

    import perseus
    r = perseus.rips(points, maxdim=2)              # Q barcode + every bad prime
    r["barcode"][1]                                 # (m, 2) array of H1 intervals, death = inf if essential
    r["certificate"]["proof"]                       # True: the barcode and the bad primes are proved
    r["certificate"]["bad_primes"]                  # {p: {"only_over_Fp": ..., "only_over_Q": ..., ...}}
    perseus.rips(points, maxdim=2, coeff=2)         # just the F_2 barcode
    perseus.rips(points, maxdim=2, mode="fast")     # random primes: faster on torsion-heavy inputs, NOT a proof
    perseus.cubical(image), perseus.flag([(0, 1, 0.5), ...]), perseus.simplices([((0, 1, 2), 3.0), ...]), perseus.cells([...])

Integral and simplicial-set inputs (see docs/TECHNIQUES.md):
    perseus.cells(cells, integral=True, image=(0, 1))     # + r["integral"]: H_*(K; Z), image of H(K_0) -> H(K_1)
    perseus.simplicial_set([(0, 0, []), (1, 0, [0, 0]), (2, 0, [1, (0, [0]), 1])], integral=True)   # RP^2
    perseus.nerve(src, tgt, identities, composition, weights), perseus.bar_construction(mult)       # B G
    perseus.magnitude(edges, lmax=5)                      # magnitude homology over Z with torsion
    perseus.pc_group(**perseus.read_pcgroup("g.txt"), maxdim=8, ring=2)   # H_*(G; Z), H^*(G; F_2) ring data
    perseus.laurent([(0, 0, []), (1, 0, [(0, "t - 1")])])   # Z[t, t^-1] coefficients: Q(t) barcode, Q[t^+-1] homology

Certificates: mode="proved" (the default; "auto" is the same) runs the Q reduction in exact rational
arithmetic, int64 where the numbers fit and big rationals where they do not, so the barcode and the bad
primes are proved (certificate "proof": True). mode="fast" uses random 62-bit primes and reports a bound on
the error probability ("log2_error_bound"); it is not a proof, and it says so: "proof": False, a "warning"
entry, and a Python warning. An entry of "bad_primes" with "every_prime_factor_of" set applies to every
prime factor of its key. See docs/CERTIFICATE.md.
"""
from ._perseus import rips, flag, cells, simplices, cubical, simplicial_set, nerve, bar_construction, magnitude, laurent, pc_group

__all__ = ["rips", "flag", "cells", "simplices", "cubical", "simplicial_set", "nerve", "bar_construction", "magnitude", "laurent", "pc_group", "read_pcgroup", "bad_primes"]
__version__ = "6.0.0"


def bad_primes(result):
    """Sorted list of the bad primes (and numbers all of whose prime factors are bad) of a Q result."""
    return sorted(result.get("certificate", {}).get("bad_primes", {}).keys())


def read_pcgroup(path):
    """Read a file written by tools/export_pcgroup.g; returns the keyword arguments of pc_group:
    mult, pcgs, relative_orders, identity."""
    tok = []
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].replace(",", " ")
            tok += [int(t) for t in line.split()]
    n, k = tok[0], tok[1]
    pcgs = tok[2:2 + k]; rel = tok[2 + k:2 + 2 * k]; e = tok[2 + 2 * k]
    flat = tok[3 + 2 * k:]
    if len(flat) != n * n:
        raise ValueError("read_pcgroup: expected an %d x %d multiplication table" % (n, n))
    return {"mult": [flat[i * n:(i + 1) * n] for i in range(n)], "pcgs": pcgs, "relative_orders": rel, "identity": e}
