#!/usr/bin/env python3
"""Fuzz the cubical engines through the Python module.

  * against the explicit complex (ground truth; small grids, absent cells allowed): F_2 and F_3 barcodes,
    the Q barcode and its bad primes, for the default (gradient) engine and the classic engine;
  * against CubicalRipser (tcripser, if installed; finite grids): positive-length F_2 intervals.

Value distributions stress ties: few integer levels, smoothed fields quantised to few levels, step
functions, +-0.0, huge and tiny magnitudes, plus degenerate shapes (axes of length 1 and 2, elongated
grids). Usage: tests/fuzz_cubical.py [trials] [seed]
"""
import sys, random
from collections import Counter
import numpy as np
import perseus

trials = int(sys.argv[1]) if len(sys.argv) > 1 else 500
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 20261118
rnd = random.Random(seed); rng = np.random.default_rng(seed)
try:
    import tcripser
except ImportError:
    tcripser = None

def bars(r):
    return Counter((k, float(b), float(e)) for k, arr in enumerate(r["barcode"]) for b, e in np.asarray(arr).reshape(-1, 2) if e > b)

def cripser_bars(A):
    out = tcripser.computePH(np.ascontiguousarray(A, dtype=float), maxdim=A.ndim, location="none")
    return Counter((int(x[0]), float(x[1]), float("inf") if x[2] >= 1.7e308 else float(x[2])) for x in out if x[2] > x[1])

def shape():
    d = rnd.choice([2, 2, 3, 3, 4])
    kind = rnd.random()
    if kind < 0.15:                                    # an axis of length 1 or 2
        s = [rnd.randint(1, {2: 40, 3: 12, 4: 6}[d]) for _ in range(d)]; s[rnd.randrange(d)] = rnd.choice([1, 2])
    elif kind < 0.25:                                  # elongated
        s = [rnd.choice([1, 2, 3]) for _ in range(d)]; s[rnd.randrange(d)] = rnd.randint(20, {2: 400, 3: 120, 4: 60}[d])
    elif kind < 0.30:
        s = [2] * d
    else:
        s = [rnd.randint(2, {2: 60, 3: 16, 4: 8}[d]) for _ in range(d)]
    return tuple(s)

def values(s):
    k = rnd.randrange(9)
    if k == 0: A = rng.random(s); name = "uniform"
    elif k == 1: A = rng.integers(0, rnd.choice([2, 3, 4]), size=s).astype(float); name = "2-4 levels"
    elif k == 2: A = rng.integers(0, rnd.choice([5, 8, 16]), size=s).astype(float); name = "5-16 levels"
    elif k == 3:
        from scipy.ndimage import gaussian_filter
        B = gaussian_filter(rng.random(s), rnd.choice([0.7, 1.5, 3])); q = rnd.choice([2, 3, 5, 10])
        A = np.floor(q * (B - B.min()) / (np.ptp(B) + 1e-12)); name = f"smooth, {q} levels"
    elif k == 4:                                       # step function along a random axis
        ax = rnd.randrange(len(s)); idx = np.indices(s)[ax]; A = (idx // max(1, s[ax] // 3)).astype(float); name = "steps"
    elif k == 5: A = rng.choice([-0.0, 0.0, 1.0, -1.0], size=s); name = "+-0"
    elif k == 6: A = rng.choice([1e-300, 1e300, -1e300, 5.0, 5.0 + 2**-40], size=s); name = "extreme magnitudes"
    elif k == 7: A = np.full(s, 7.0); A[tuple(rng.integers(0, n) for n in s)] = 1.0; name = "constant with one dip"
    else: A = rng.normal(size=s).round(1); name = "rounded normal"
    return A, name

res = Counter(); fails = []
def check(ok, what):
    res["pass" if ok else "fail"] += 1
    if not ok: fails.append(what); print("FAIL", what, flush=True)

for s in [(0, 5), (3, 0, 2), (2, 2, 0, 2)]:          # axes of length 0: the empty complex
    for e in ("auto", "classic", "morse", "explicit"):
        check(bars(perseus.cubical(np.zeros(s), engine=e)) == Counter(), f"empty {s} engine {e}")

for t in range(trials):
    s = shape(); A, name = values(s)
    tag = f"#{t} {s} {name}"
    small = A.size <= 4000
    if small and rnd.random() < 0.4:                   # absent cells: explicit complex only
        B = A.copy(); B[rng.random(s) < rnd.choice([0.1, 0.3, 0.6, 1.0])] = np.nan
        for p in (2, 3):
            check(bars(perseus.cubical(B, coeff=p)) == bars(perseus.cubical(B, coeff=p, engine="explicit")), f"{tag} absent F_{p} default == explicit")
        a, b = perseus.cubical(B), perseus.cubical(B, engine="explicit")
        check(bars(a) == bars(b) and perseus.bad_primes(a) == perseus.bad_primes(b), f"{tag} absent Q default == explicit")
        continue
    if small:
        for p in (2, 3):
            e = bars(perseus.cubical(A, coeff=p, engine="explicit"))
            check(bars(perseus.cubical(A, coeff=p)) == e, f"{tag} F_{p} default == explicit")
            check(bars(perseus.cubical(A, coeff=p, engine="classic")) == e, f"{tag} F_{p} classic == explicit")
        a, b = perseus.cubical(A), perseus.cubical(A, engine="explicit")
        check(bars(a) == bars(b) and perseus.bad_primes(a) == perseus.bad_primes(b), f"{tag} Q default == explicit")
    m = bars(perseus.cubical(A, coeff=2))
    if tcripser is not None:
        check(m == cripser_bars(A), f"{tag} F_2 default == CubicalRipser")
    if not small:
        check(m == bars(perseus.cubical(A, coeff=2, engine="classic")), f"{tag} F_2 default == classic")

print(f"\n{res['pass']} passed, {res['fail']} failed (seed {seed}, {trials} trials, CubicalRipser {'on' if tcripser else 'not installed'})")
for f in fails[:30]: print("  ", f)
sys.exit(1 if res["fail"] else 0)
