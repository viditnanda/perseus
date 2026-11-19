#!/usr/bin/env python3
"""Independent check of `perseus --torsion`: torsion of H_{q}(K_t; Z) at every filtration value t, by FLINT
Smith forms of the boundary matrix d_{q+1} restricted to the cells born by t.

usage: tests/torsion_flint.py kind n m seed [d]   (needs python-flint)
  kind = cells: random filtered 2-level chain complex over Z (n 1-cells, m 2-cells; tied values if seed is even)
  kind = lm:  d-dimensional Linial-Meshulam process (complete (d-1)-skeleton at value 1, then m random d-cells)
  kind = mix: 2-complex with distinct interleaved values (edges random in [0,1], each triangle after its edges)
"""
import os, tempfile, sys, itertools, json, subprocess, random
import flint

kind, n, m, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
d = int(sys.argv[5]) if len(sys.argv) > 5 else 2
rng = random.Random(seed)
lines = []
if kind == "lm":
    faces = list(itertools.combinations(range(n), d))
    tops = list(itertools.combinations(range(n), d + 1)); rng.shuffle(tops); tops = tops[:m]
    for f in faces: lines.append((f, 1.0))
    for i, t in enumerate(tops): lines.append((t, float(i + 2)))
elif kind == "cells":
    d = 2
    # n 1-cells (boundary 0), m 2-cells with sparse random integer boundaries; tied rows if seed even
    tied = seed % 2 == 0
    b1 = [0.0 if tied else float(f"{rng.random():.12f}") for _ in range(n)]
    coefs = [1, -1, 2, -2, 3, 4, -4, 6, 8, 9, -9, 12, 27]
    cl = []
    for i in range(n): cl.append(f"1 {b1[i]:.12f}")
    cellsbd = []
    for j in range(m):
        k = rng.choice([1, 2, 2, 3])
        fs = rng.sample(range(n), k); bd = [(f, rng.choice(coefs)) for f in fs]
        v = (1.0 + j) if tied else max(b1[f] for f, _ in bd) + rng.random()
        cellsbd.append((bd, float(f"{v:.12f}")))
    path = os.path.join(tempfile.gettempdir(), f"chk_cells_{n}_{m}_{seed}.txt")
    with open(path, "w") as f:
        for l in cl: f.write(l + "\n")
        for bd, v in cellsbd: f.write(f"2 {v:.12f} " + " ".join(f"{a} {c}" for a, c in bd) + "\n")
    out = subprocess.run([os.environ.get("PERSEUS", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "perseus")), "cells", path, "--maxdim", "2", "--torsion", "--print", "json"], capture_output=True, text=True)
    J = json.loads(out.stdout); T = J["torsion"]; assert T["complete"], T
    bad = sorted(int(b["p"]) for b in J["bad_primes"])
    steps = next((x["steps"] for x in T["degrees"] if x["degree"] == 1), [])
    def at(t):
        cur = []
        for s_ in steps:
            if s_["value"] <= t: cur = [int(x) for x in s_["invariant_factors"]]
        return cur
    nm = 0; nchk = 0
    for t in sorted({v for _, v in cellsbd}):
        C = [bd for bd, v in cellsbd if v <= t]
        M = flint.fmpz_mat(n, len(C))
        for j, bd in enumerate(C):
            for a, c in bd: M[a, j] = c
        S = M.snf(); diag = [abs(int(S[i, i])) for i in range(min(S.nrows(), S.ncols()))]
        tors = sorted(x for x in diag if x > 1); nchk += 1
        if tors != at(t):
            nm += 1
            if nm <= 5: print(f"  MISMATCH at t={t}: smith {tors} vs perseus {at(t)}")
    mx = max((max([int(f) for f in s_["invariant_factors"]] or [1]) for s_ in steps), default=1)
    print(f"cells n={n} m={m} seed={seed} ({'tied' if tied else 'distinct'}): {nchk} groups checked, {nm} mismatches; {len(steps)} steps; largest invariant factor {mx}; bad primes {bad}")
    sys.exit(1 if nm else 0)
else:
    d = 2
    E = {e: rng.random() for e in itertools.combinations(range(n), 2)}
    E = {e: float(f"{w:.12f}") for e, w in E.items()}
    for e, w in E.items(): lines.append((e, w))
    tris = list(itertools.combinations(range(n), 3)); rng.shuffle(tris)
    for t in tris[:m]:
        me = max(E[(t[0], t[1])], E[(t[0], t[2])], E[(t[1], t[2])])
        lines.append((t, float(f"{max(rng.random() * 1.3, me + 1e-7 * (1 + rng.random())):.12f}")))
path = os.path.join(tempfile.gettempdir(), f"chk_{kind}{d}_{n}_{m}_{seed}.txt")
with open(path, "w") as f:
    for s, v in lines: f.write(" ".join(map(str, s)) + f" {v:.12f}\n")

# value of every simplex (faces: min of cofaces), as Perseus reads it
val = {}
for s, v in lines:
    for k in range(1, len(s) + 1):
        for f in itertools.combinations(s, k): val[f] = min(val.get(f, float("inf")), v)

out = subprocess.run([os.environ.get("PERSEUS", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "perseus")), "simplices", path, "--maxdim", str(d), "--torsion", "--print", "json"],
                     capture_output=True, text=True)
J = json.loads(out.stdout)
T = J["torsion"]
assert T["complete"], T
bad = sorted(int(b["p"]) for b in J["bad_primes"])

def predicted(q):
    """value -> invariant factors (list of ints) from the Perseus steps (piecewise constant)"""
    steps = next((x["steps"] for x in T["degrees"] if x["degree"] == q), [])
    def at(t):
        cur = []
        for s in steps:
            if s["value"] <= t: cur = [int(x) for x in s["invariant_factors"]]
        return cur
    return at

nmis = 0; nchk = 0; primes_seen = set()
for q in range(0, d):
    rows = sorted([s for s in val if len(s) == q + 1], key=lambda s: (val[s], s))
    cols = sorted([s for s in val if len(s) == q + 2], key=lambda s: (val[s], s))
    if not cols: continue
    ri = {s: i for i, s in enumerate(rows)}
    at = predicted(q)
    values = sorted({val[c] for c in cols})
    for t in values:
        C = [c for c in cols if val[c] <= t]
        M = flint.fmpz_mat(len(rows), len(C))
        for j, c in enumerate(C):
            for i in range(len(c)):
                M[ri[c[:i] + c[i + 1:]], j] = (-1) ** i
        S = M.snf()
        diag = [abs(int(S[i, i])) for i in range(min(S.nrows(), S.ncols()))]
        tors = sorted(x for x in diag if x > 1)
        for x in tors:
            for p, e in flint.fmpz(x).factor(): primes_seen.add(int(p))
        nchk += 1
        if tors != at(t):
            nmis += 1
            if nmis <= 5: print(f"  MISMATCH H_{q} at t={t}: smith {tors} vs perseus {at(t)}")
steps = sum(len(x["steps"]) for x in T["degrees"])
maxtors = max((max([int(f) for f in s["invariant_factors"]] or [1]) for x in T["degrees"] for s in x["steps"]), default=1)
print(f"{kind} d={d} n={n} m={m} seed={seed}: {nchk} groups checked, {nmis} mismatches; {steps} steps; "
      f"largest invariant factor {maxtors}; bad primes {bad}; primes in absolute torsion {sorted(primes_seen)}")
sys.exit(1 if nmis else 0)
