#!/usr/bin/env python3
"""Independent check: barcodes over Q and F_p of a 'simplices' file (faces get the min of coface values),
by plain column reduction in filtration order (value, dim, lexicographic). Value-indexed barcodes
(zero-length bars dropped). Prints which primes change the barcode."""
import sys, itertools
from fractions import Fraction

lines = [l.split() for l in open(sys.argv[1])]
primes = [int(x) for x in sys.argv[2].split(",")]
val = {}
for t in lines:
    s = tuple(sorted(map(int, t[:-1]))); v = float(t[-1])
    for k in range(1, len(s) + 1):
        for f in itertools.combinations(s, k):
            val[f] = min(val.get(f, float("inf")), v)
cells = sorted(val, key=lambda f: (val[f], len(f), f))
pos = {c: i for i, c in enumerate(cells)}

def bnd(c):
    if len(c) == 1: return {}
    return {pos[c[:i] + c[i + 1:]]: (-1) ** i for i in range(len(c))}

def reduce(field):
    one = field(1); low_of = {}; pairs = []
    for j, c in enumerate(cells):
        col = {r: field(x) for r, x in bnd(c).items()}
        col = {r: x for r, x in col.items() if x != 0}
        while col:
            l = max(col)
            if l not in low_of: break
            k, kc = low_of[l]
            f = col[l] / kc[l]
            for r, x in kc.items():
                y = col.get(r, 0) - f * x
                if y == 0: col.pop(r, None)
                else: col[r] = y
        if col:
            l = max(col); low_of[l] = (j, col); pairs.append((l, j))
    paired = {l for l, _ in pairs} | {j for _, j in pairs}
    bars = sorted((len(cells[l]) - 1, val[cells[l]], val[cells[j]]) for l, j in pairs if val[cells[l]] < val[cells[j]])
    bars += sorted((len(cells[i]) - 1, val[cells[i]], float("inf")) for i in range(len(cells)) if i not in paired)
    return sorted(bars)

class Fp:
    p = 2
    def __init__(s, x): s.x = x % Fp.p
    def __sub__(s, o): return Fp(s.x - o.x)
    def __rsub__(s, o): return Fp(o - s.x)
    def __mul__(s, o): return Fp(s.x * o.x)
    def __rmul__(s, o): return Fp(s.x * o.x)
    def __truediv__(s, o): return Fp(s.x * pow(o.x, -1, Fp.p))
    def __eq__(s, o): return s.x == (o.x if isinstance(o, Fp) else o % Fp.p)
    def __ne__(s, o): return not s == o

q = reduce(Fraction)
for p in primes:
    Fp.p = p
    b = reduce(Fp)
    print(p, "bad" if b != q else "harmless", "(%d bars differ)" % len(set(b) ^ set(q)))
