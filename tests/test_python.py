#!/usr/bin/env python3
"""The Python module against the command line on the same inputs: barcodes, bad primes and certificate
fields, in every certificate mode. Needs the module installed (cd python && pip install .) and ./perseus."""
import itertools, random, json, os, subprocess, sys, tempfile
from collections import Counter
import numpy as np
import perseus

HERE = os.path.dirname(os.path.abspath(__file__))
PERSEUS = os.environ.get("PERSEUS", os.path.join(HERE, "..", "perseus"))
FIX = os.path.join(HERE, "fixtures")
TMP = tempfile.mkdtemp(prefix="p6py_")
results = Counter(); failures = []

def check(name, ok, detail=""):
    results["pass" if ok else "fail"] += 1
    if not ok:
        failures.append(name + (": " + str(detail)[:300] if detail else ""))
        print("FAIL", name, str(detail)[:300], flush=True)

def key(b, d): return ("%.9g" % b, "inf" if d is None or d == float("inf") else "%.9g" % d)

def py_bars(r):
    out = Counter()
    for k, arr in enumerate(r["barcode"]):
        for b, d in arr: out[(k,) + key(b, d)] += 1
    return out

def cli(args):
    p = subprocess.run([PERSEUS] + args + ["--print", "json"], capture_output=True, text=True, timeout=600)
    if p.returncode: raise RuntimeError(p.stderr)
    j = json.loads(p.stdout); out = Counter()
    for k, lst in enumerate(j["barcode"]):
        for b, d in lst: out[(k,) + key(b, d)] += 1
    return out, j

def compare(name, r, args, mode="auto"):
    a = py_bars(r); b, j = cli(args + (["--mode", mode] if "--coeff" not in args else []))
    check(f"{name} barcode", a == b, f"{dict(a - b)} {dict(b - a)}")
    if "certificate" in r:
        c = r["certificate"]
        check(f"{name} bad primes", sorted(c["bad_primes"]) == sorted(int(x["p"]) for x in j["bad_primes"]), f"{sorted(c['bad_primes'])} vs {[x['p'] for x in j['bad_primes']]}")
        check(f"{name} mode", c["mode"] == j["mode"], f"{c['mode']} vs {j['mode']}")
        check(f"{name} complete", c["complete"] == j["complete"])

rng = np.random.default_rng(3)

# rips: points (implicit and sparse routes), square and condensed distances
for t in range(6):
    n = int(rng.integers(20, 120)); X = rng.random((n, int(rng.integers(2, 4)))); md = int(rng.integers(1, 3))
    path = f"{TMP}/pts{t}.txt"; np.savetxt(path, X)
    for mode in ("auto", "fast", "proved"):
        r = perseus.rips(points=X, maxdim=md, threshold=0.5, mode=mode, seed=7)
        compare(f"rips points n={n} H<={md} {mode}", r, ["points", path, "--maxdim", str(md), "--threshold", "0.5", "--seed", "7"], mode)
    for eng in ("implicit", "sparse"):
        r = perseus.rips(points=X, maxdim=md, threshold=0.5, coeff=2, engine=eng)
        compare(f"rips points {eng} F2", r, ["points", path, "--maxdim", str(md), "--threshold", "0.5", "--coeff", "2", "--engine", eng])
    D = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
    a = py_bars(perseus.rips(distances=D, maxdim=md, threshold=0.5)); b = py_bars(perseus.rips(distances=D[np.tril_indices(n, -1)], maxdim=md, threshold=0.5))
    c = py_bars(perseus.rips(points=X, maxdim=md, threshold=0.5))
    check(f"rips square == condensed == points (n={n})", a == b == c)

# a sparse cloud large enough for the automatic sparse route
X = rng.random((3000, 3)); np.savetxt(f"{TMP}/big.txt", X)
compare("rips 3000 points auto route", perseus.rips(points=X, maxdim=2, threshold=0.06), ["points", f"{TMP}/big.txt", "--maxdim", "2", "--threshold", "0.06"])

# flag, with and without vertex births
for t in range(6):
    n = int(rng.integers(5, 14)); vb = [float(rng.choice([0, 0, 1, 2])) for _ in range(n)]
    E = [(i, j, max(vb[i], vb[j]) + float(rng.integers(0, 6))) for i in range(n) for j in range(i + 1, n) if rng.random() < 0.6]
    for births in (vb, [0.0] * n):
        path = f"{TMP}/f{t}.txt"
        with open(path, "w") as fh:
            for i in range(n): fh.write(f"{i} {births[i]}\n")
            for i, j, w in E: fh.write(f"{i} {j} {max(w, births[i], births[j])}\n")
        EE = [(i, j, max(w, births[i], births[j])) for i, j, w in E]
        r = perseus.flag(EE, vertex_births=births, maxdim=2)
        compare(f"flag n={n} births={'yes' if any(births) else 'no'}", r, ["flag", path, "--maxdim", "2"])

# cells: the certificate fixtures, including a bad prime above 2^64 and a group entry
def read_cells(path):
    cells = []
    for line in open(path):
        tok = line.split()
        if not tok: continue
        d, b = int(tok[0]), float(tok[1]); rest = list(map(int, tok[2:]))
        cells.append((d, b, list(zip(rest[::2], rest[1::2]))))
    return cells
for f in ("crt_alias_0.cells", "crt_alias_1.cells", "crt_group.cells", "crt_group_split.cells", "crt_group_degree.cells", "crt_shared_degree.cells"):
    cells = read_cells(os.path.join(FIX, f))
    for mode in ("auto", "fast", "proved"):
        r = perseus.cells(cells, mode=mode, seed=1)
        compare(f"cells {f} {mode}", r, ["cells", os.path.join(FIX, f), "--seed", "1"], mode)
        r2 = perseus.cells(cells, mode=mode)
        _, j = cli(["cells", os.path.join(FIX, f), "--mode", mode])
        mine = sorted((p, d["every_prime_factor_of"], d["first_degree"]) for p, d in r2["certificate"]["bad_primes"].items())
        theirs = sorted((int(x["p"]), bool(x.get("every_prime_factor_of")), x["first_degree"]) for x in j["bad_primes"])
        check(f"cells {f} {mode} keys, groups and degrees as on the command line", mine == theirs and 0 not in r2["certificate"]["bad_primes"], f"{mine} vs {theirs}")
# ---- the default certificate is a proof, also where int64 overflows ------------------------------------------
# A Linial-Meshulam process near the threshold (all edges of K_60, then 1700 random triangles): int64 overflows,
# and the default must still prove the barcode and the bad primes (exact big rationals), not fall back to random
# primes. Fast mode must agree, and must say it is not a proof.
_rng = np.random.default_rng(5)
_tris = list(itertools.combinations(range(60), 3)); _idx = _rng.permutation(len(_tris))[:1700]
_lm = [((a, b), 1.0) for a, b in itertools.combinations(range(60), 2)] + [(_tris[i], float(t + 2)) for t, i in enumerate(_idx)]
r = perseus.simplices(_lm)
c = r["certificate"]
check("default on a torsion-heavy complex: proved, a proof", c["mode"] == "proved" and c["proof"] and c["complete"] and "big rationals" in c["arithmetic"], (c["mode"], c.get("arithmetic")))
check("default on a torsion-heavy complex: bad primes 2, 3, 922457715037", sorted(c["bad_primes"]) == [2, 3, 922457715037], sorted(c["bad_primes"]))
r_auto = perseus.simplices(_lm, mode="auto")
check("mode='auto' is the same as proved", r_auto["certificate"]["mode"] == "proved" and py_bars(r_auto) == py_bars(r) and sorted(r_auto["certificate"]["bad_primes"]) == sorted(c["bad_primes"]))
import warnings
with warnings.catch_warnings(record=True) as _w:
    warnings.simplefilter("always")
    rf = perseus.simplices(_lm, mode="fast")
cf = rf["certificate"]
check("fast mode: same barcode and bad primes, labelled as no proof",
      py_bars(rf) == py_bars(r) and sorted(cf["bad_primes"]) == sorted(c["bad_primes"]) and cf["mode"] == "fast" and not cf["proof"]
      and "not a proof" in cf["warning"] and any("not a proof" in str(x.message) for x in _w), (cf["mode"], cf.get("warning"), [str(x.message) for x in _w]))
_path = os.path.join(TMP, "lm60.txt")
with open(_path, "w") as _fh:
    for s, b in _lm: _fh.write(" ".join(map(str, s)) + f" {int(b)}\n")
_p = subprocess.run([PERSEUS, "simplices", _path, "--print", "json"], capture_output=True, text=True, timeout=600); _j = json.loads(_p.stdout)
check("command line default: proof true, no warning", _j["mode"] == "proved" and _j["proof"] is True and "warning" not in _j and "WARNING" not in _p.stderr)
_p = subprocess.run([PERSEUS, "simplices", _path, "--print", "json", "--mode", "fast"], capture_output=True, text=True, timeout=600); _j = json.loads(_p.stdout)
check("command line fast: proof false, warning in the JSON and on stderr", _j["mode"] == "fast" and _j["proof"] is False and "not a proof" in _j["warning"] and "WARNING: --mode fast" in _p.stderr)

r = perseus.cells(read_cells(os.path.join(FIX, "crt_alias_1.cells")))
big = [p for p, d in r["certificate"]["bad_primes"].items() if d["probable_prime"]]
check("cells: 98-bit bad prime reported", big == [294719631262177595678027596519], big)
r = perseus.cells(read_cells(os.path.join(FIX, "crt_group.cells")))
g = [p for p, d in r["certificate"]["bad_primes"].items() if d["every_prime_factor_of"]]
check("cells: group entry", g == [36028797018976327 * 1152921504606848017], g)
r = perseus.cells(read_cells(os.path.join(FIX, "crt_alias_1.cells")), mode="fast", seed=3)
c = r["certificate"]
check("fast mode with a supplied seed: bound and seed reported, not certified", c["log2_error_bound"] < -64 and int(c["seed"], 16) == 3
      and c["seed_supplied"] and not c["barcode_certified"] and not c["complete"])
r = perseus.cells(read_cells(os.path.join(FIX, "crt_alias_1.cells")), mode="fast"); c = r["certificate"]
check("fast mode, fresh seed: certified", c["barcode_certified"] and c["complete"] and len(c["seed"]) == 64 and not c["seed_supplied"])
r2 = perseus.cells(read_cells(os.path.join(FIX, "crt_alias_1.cells")), mode="fast", seed=c["seed"])
check("fast mode: the reported seed reproduces the run", r2["certificate"]["moduli"] == c["moduli"])
r = perseus.cells([(0, 0.0, []), (1, 1.0, []), (2, 2.0, [(1, 12)])])
check("Moore space M(Z/12,1): bad primes [2, 3]", perseus.bad_primes(r) == [2, 3], perseus.bad_primes(r))

# simplices: RP^2
T = [(0,1,2),(0,2,3),(0,3,4),(0,4,5),(0,5,1),(1,2,4),(2,3,5),(3,4,1),(4,5,2),(5,1,3)]
r = perseus.simplices([(list(t), 0.0) for t in T])
check("RP2: bad prime 2", perseus.bad_primes(r) == [2], perseus.bad_primes(r))
compare("simplices RP2", r, ["simplices", os.path.join(HERE, "..", "examples", "rp2_simplices.txt")])

# cubical: implicit engine vs the command line, periodic axes and absent cells
for t in range(10):
    d = int(rng.integers(1, 4)); shape = tuple(int(rng.integers(1, 7)) for _ in range(d))
    A = rng.integers(0, 6, size=shape).astype(float)
    if rng.random() < 0.4: A[rng.random(shape) < 0.15] = np.nan
    per = [bool(rng.random() < 0.35) for _ in range(d)]
    path = f"{TMP}/c{t}.npy"; np.save(path, A)
    pa = ["--periodic", ",".join(str(i) for i in range(d) if per[i])] if any(per) else []
    for coeff in ("Q", 3):
        r = perseus.cubical(A, periodic=per, coeff=coeff)
        compare(f"cubical {shape} per={per} {coeff}", r, ["npy", path, "--coeff", str(coeff)] + pa)
    check(f"cubical {shape} implicit == explicit", py_bars(perseus.cubical(A, periodic=per)) == py_bars(perseus.cubical(A, periodic=per, engine="explicit")))
    check(f"cubical {shape} classic engine == default", py_bars(perseus.cubical(A, periodic=per, engine="classic")) == py_bars(perseus.cubical(A, periodic=per)))
    if 2 <= d <= 4 and not any(per):
        check(f"cubical {shape} morse engine == explicit", py_bars(perseus.cubical(A, engine="morse")) == py_bars(perseus.cubical(A, engine="explicit")))

# extra fields
r = perseus.rips(points=rng.random((40, 2)), maxdim=1, coeff=2, fp=[3, 5])
check("fp fields", sorted(r["fields"]) == [3, 5])

# ---- integral outputs, simplicial sets, nerves, magnitude homology ------------------------------------
def groups(h): return [(g["free"], [int(t) for t in g["torsion"]]) for g in h]
rp2 = [(0, 0, []), (1, 0, []), (2, 1, [(1, 2)])]            # a circle at 0, a degree-2 disc at 1
r = perseus.cells(rp2, integral=True, image=(0, 1))
check("cells RP2: integral homology Z, Z/2, 0", groups(r["integral"]["homology"]) == [(1, []), (0, [2]), (0, [])], r["integral"])
check("cells RP2: image of H_1(K_0) -> H_1(K_1) is Z/2", groups(r["integral"]["image"]["groups"])[1] == (0, [2]))
check("cells RP2: F_2 barcode of the bad prime", sorted(r["integral"]["bad_prime_barcodes"]) == [2] and len(r["integral"]["bad_prime_barcodes"][2][2]) == 1)
try: perseus.cells([(0, 0, []), (1, 0, [(0, 1)]), (2, 0, [(1, 1)])]); check("cells: d o d != 0 rejected", False)
except ValueError: check("cells: d o d != 0 rejected", True)
# the RP2 of 6 vertices with a projective line born first: the image of H_1 is Z/2
T = [(0,1,2),(0,2,3),(0,3,4),(0,4,5),(0,5,1),(1,2,4),(2,3,5),(3,4,1),(4,5,2),(5,1,3)]
sx = {}
for t in T:
    for k in range(1, 4):
        for f in itertools.combinations(sorted(t), k): sx[f] = 1.0
for e in [(0, 1), (1, 3), (0, 3)]: sx[e] = 0.0
for v in range(6): sx[(v,)] = 0.0
order = sorted(sx, key=lambda f: (sx[f], len(f), f)); ix = {f: i for i, f in enumerate(order)}
cl = [(len(f) - 1, sx[f], [(ix[f[:i] + f[i + 1:]], (-1) ** i) for i in range(len(f))] if len(f) > 1 else []) for f in order]
r = perseus.cells(cl, integral=True, image=(0, 1))
check("RP2 (6 vertices): reduction 31 -> at most 9 cells", r["reduction"]["cells"] == 31 and r["reduction"]["after"] <= 9, r["reduction"])
check("RP2 (6 vertices): image of H_1(projective line) is Z/2", groups(r["integral"]["image"]["groups"])[1] == (0, [2]), r["integral"]["image"])
r0 = perseus.cells(cl, reduce=False, zero_length=True); r1 = perseus.cells(cl, zero_length=True)
check("RP2 (6 vertices): reduction keeps the barcode with zero-length bars", py_bars(r0) == py_bars(r1))
for f in ("crt_alias_0.cells", "crt_group.cells", "crt_unit_0.cells", "crt_unit_1.cells"):
    cells = read_cells(os.path.join(FIX, f))
    a, b = perseus.cells(cells, reduce=False, zero_length=True), perseus.cells(cells, zero_length=True)
    check(f"cells {f}: reduction keeps the barcode and bad primes", py_bars(a) == py_bars(b) and perseus.bad_primes(a) == perseus.bad_primes(b))

r = perseus.simplicial_set([(0, 0, []), (1, 0, [0, 0]), (2, 1, [1, (0, [0]), 1])], integral=True, image=(0, 1))
check("simplicial set: minimal RP2", groups(r["integral"]["homology"]) == [(1, []), (0, [2]), (0, [])] and groups(r["integral"]["image"]["groups"])[1] == (0, [2]))
r = perseus.simplicial_set([(0, 0, []), (1, 0, [0, 0]), (1, 0, [0, 0]), (1, 0, [0, 0]), (2, 0, [1, 3, 2]), (2, 0, [2, 3, 1])], integral=True)
check("simplicial set: minimal torus", groups(r["integral"]["homology"]) == [(1, []), (2, []), (1, [])], r["integral"])
try: perseus.simplicial_set([(0, 0, []), (0, 0, []), (1, 0, [0, 1]), (2, 0, [2, 2, 2])]); check("simplicial set: bad identities rejected", False)
except ValueError: check("simplicial set: bad identities rejected", True)
for n, expect in ((2, [(1, []), (0, [2]), (0, []), (0, [2]), (0, []), (0, [2])]), (3, [(1, []), (0, [3]), (0, []), (0, [3]), (0, []), (0, [3])])):
    r = perseus.bar_construction([[(a + b) % n for b in range(n)] for a in range(n)], 0, maxdim=5)
    check(f"B Z/{n}: H_k up to 5", groups(r["integral"]["homology"]) == expect, groups(r["integral"]["homology"]))
perms = list(itertools.permutations(range(3))); pid = {q: i for i, q in enumerate(perms)}
mult = [[pid[tuple(a[b[i]] for i in range(3))] for b in perms] for a in perms]
r = perseus.bar_construction(mult, pid[(0, 1, 2)], maxdim=4)
check("B S_3: Z, Z/2, 0, Z/6, 0", groups(r["integral"]["homology"]) == [(1, []), (0, [2]), (0, []), (0, [6]), (0, [])], groups(r["integral"]["homology"]))
# the nerve of a weighted poset is its weighted order complex
for t in range(6):
    n = 6; rel = [(i, j) for i in range(n) for j in range(i + 1, n) if rng.random() < 0.45]
    less = {(i, j) for i, j in rel}
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if (i, k) in less and (k, j) in less: less.add((i, j))
    wv = [float(rng.integers(0, 3)) for _ in range(n)]
    mor = [(i, i) for i in range(n)] + sorted(less); mi = {m: q for q, m in enumerate(mor)}
    w = [wv[i] for i in range(n)] + [float(max(wv[i], wv[j]) + rng.integers(0, 2)) for i, j in sorted(less)]
    comp = {(mi[(b, c)], mi[(a, b)]): mi[(a, c)] for (a, b) in mor for (b2, c) in mor if b == b2}
    r = perseus.nerve([m[0] for m in mor], [m[1] for m in mor], list(range(n)), comp, weights=w, maxdim=2)
    # order complex: a chain is born at the largest weight of a relation between its elements (and its vertices)
    simp = []
    for k in range(1, 5):
        for ch in itertools.combinations(range(n), k):
            if all((ch[i], ch[i + 1]) in less for i in range(k - 1)):
                simp.append((list(ch), max([wv[x] for x in ch] + [w[mi[(ch[i], ch[j])]] for i in range(k) for j in range(i + 1, k)])))
    s2 = perseus.simplices(simp, maxdim=2)
    check(f"nerve of a weighted poset == order complex ({t})", py_bars(r) == py_bars(s2), f"{py_bars(r)} vs {py_bars(s2)}")
# a category with a non-identity endomorphism: the idempotent monoid {1, e} has a contractible nerve
r = perseus.bar_construction([[0, 1], [1, 1]], 0, maxdim=4)
check("B of the monoid {1, e = e^2} is contractible", groups(r["integral"]["homology"]) == [(1, [])] + [(0, [])] * 4)

r = perseus.magnitude([(i, (i + 1) % 5) for i in range(5)], lmax=5)
c5 = {k: v["free"] for k, v in r["table"].items() if v["free"]}
check("magnitude C_5 (Hepworth-Willerton)", c5 == {(0, 0): 5, (1, 1): 10, (2, 2): 10, (2, 3): 10, (3, 3): 10, (3, 4): 30, (4, 4): 10, (4, 5): 50, (5, 5): 10}, c5)
r = perseus.magnitude([(0, 1), (1, 2), (1, 3), (3, 4)], lmax=4)
check("magnitude of a tree is diagonal and torsion-free", all(k == l and not v["torsion"] for (k, l), v in r["table"].items() if v["free"] or v["torsion"]))
r = perseus.magnitude([(i, (i + 1) % 4) for i in range(4)], directed=True, lmax=4)
check("magnitude: the directed 4-cycle", all(v["free"] in (0, 4) and not v["torsion"] for v in r["table"].values()))
# command line: simplicial sets and magnitude agree with the module
path = f"{TMP}/rp2.sset"; open(path, "w").write("0 0\n1 0 0 0\n2 1 1 0@0 1\n")
_, j = cli(["sset", path, "--integral", "--image", "0,1"])
check("command line sset: integral homology and image", [(g["free"], g["torsion"]) for g in j["integral_homology"]] == [(1, []), (0, [2]), (0, [])]
      and j["integral_image"]["groups"][1]["torsion"] == [2])

# ---- Laurent coefficients Z[t, t^-1] ----------------------------------------------------------------------
def lbars(B): return [sorted((float(b), float(d)) for b, d in np.asarray(x).reshape(-1, 2)) for x in B]
tref = [(0, 0, []), (1, 0, [(0, "t - 1")]), (1, 0, [(0, {1: 1, 0: -1})]), (2, 0, [(1, "t^2 - t + 1"), (2, "-t^2 + t - 1")])]
r = perseus.laurent(tref)
check("laurent trefoil: Alexander module", [(g["free"], g["torsion"]) for g in r["homology"]] == [(0, ["t - 1"]), (0, ["t^2 - t + 1"]), (0, [])], r["homology"])
check("laurent trefoil: bad irreducibles t - 1, t^2 - t + 1", sorted(b["polynomial"] for b in r["bad_irreducibles"]) == ["t - 1", "t^2 - t + 1"] and all(b["irreducible"] for b in r["bad_irreducibles"]))
check("laurent trefoil: empty Q(t) barcode, complete", all(len(np.asarray(x)) == 0 for x in r["barcode"]) and r["complete"] and not r["bad_primes"] and not r["homology_bad_primes"])
r2 = perseus.laurent([(0, 0, []), (1, 0, [(0, "1 - t^-1")]), (1, 0, [(0, "t - 1")]), (2, 0, [(1, "t^3 - t^2 + t"), (2, "-t^2 + t - 1")])])
check("laurent: units t^k do not change the answer", r2["homology"] == r["homology"], r2["homology"])
r = perseus.laurent([(0, 0, []), (1, 0, [(0, 2)]), (1, 0, [(0, "t + 1")])])
check("laurent: L/(2, t + 1) is zero over Q[t^+-1], F_2[t^+-1]/(t + 1) at p = 2",
      [(g["free"], g["torsion"]) for g in r["homology"]] == [(0, []), (1, [])] and list(r["homology_bad_primes"]) == [2]
      and r["homology_bad_primes"][2][0]["torsion"] == ["t + 1"], (r["homology"], r["homology_bad_primes"]))
try: perseus.laurent([(0, 0, []), (1, 0, [(0, "t - 1")]), (2, 0, [(1, 1)])]); check("laurent: d o d != 0 rejected", False)
except ValueError: check("laurent: d o d != 0 rejected", True)
cells_f = [(1, 0, []), (2, 1, [(0, "2t^3 - 6t^2 + 6t - 4")])]
r = perseus.laurent(cells_f)
check("laurent: filtered, bars [0,1) over Q(t)", lbars(r["barcode"]) == [[], [(0.0, 1.0)], []])
check("laurent: bad irreducibles t - 2, t^2 - t + 1 and bad prime 2", [b["polynomial"] for b in r["bad_irreducibles"]] == ["t - 2", "t^2 - t + 1"] and list(r["bad_primes"]) == [2])
path = f"{TMP}/f.lau"; open(path, "w").write("1 0\n2 1 0:2t^3-6t^2+6t-4\n")
p = subprocess.run([PERSEUS, "laurent", path, "--print", "json"], capture_output=True, text=True); j = json.loads(p.stdout)
inf = lambda B: [[[a, float("inf") if b is None else b] for a, b in x] for x in B]
check("laurent: command line agrees with the module",
      inf(j["barcode"]) == [[list(map(float, b)) for b in np.asarray(x).reshape(-1, 2)] for x in r["barcode"]]
      and [b["polynomial"] for b in j["bad_irreducibles"]] == [b["polynomial"] for b in r["bad_irreducibles"]]
      and j["homology"] == r["homology"] and sorted(map(int, j["bad_primes"])) == sorted(r["bad_primes"]), j)

# ---- finite solvable groups (pc_group) --------------------------------------------------------------------
PCF = os.path.join(FIX, "pcgroups")
r = perseus.pc_group(**perseus.read_pcgroup(os.path.join(PCF, "g8_4.txt")), maxdim=6, ring=2)
check("pc_group Q8: integral homology", [(h["free"], h["torsion"]) for h in r["homology"]] == [(1, []), (0, [2, 2]), (0, []), (0, [8]), (0, []), (0, [2, 2]), (0, [])], r["homology"])
c = r["cohomology"]
check("pc_group Q8: mod-2 dims and generators", c["dims"] == [1, 2, 2, 1, 1, 2, 2] and c["generators_by_degree"][1:5] == [2, 0, 0, 1], c)
x2, xy, y2 = c["products"][(1, 0, 1, 0)], c["products"][(1, 0, 1, 1)], c["products"][(1, 1, 1, 1)]
check("pc_group Q8: x^2 + xy + y^2 = 0 in H^2(Q8; F_2)", [(a + b + d) % 2 for a, b, d in zip(x2, xy, y2)] == [0, 0], (x2, xy, y2))
r = perseus.pc_group(**perseus.read_pcgroup(os.path.join(PCF, "g24_3.txt")), maxdim=5)
check("pc_group SL(2,3): mixed relative orders, against HAP", [h["torsion"] for h in r["homology"]][1:] == [[3], [], [24], [], [3]], r["homology"])
check("pc_group: relative orders are inferred", perseus.pc_group(**{**perseus.read_pcgroup(os.path.join(PCF, "g24_3.txt")), "relative_orders": []}, maxdim=3)["relative_orders"] == [3, 2, 2, 2])
try:
    g = perseus.read_pcgroup(os.path.join(PCF, "g8_3.txt")); perseus.pc_group(g["mult"], g["pcgs"][::-1], g["relative_orders"], g["identity"]); ok = False
except ValueError: ok = True
check("pc_group: a reversed pcgs of D8 is rejected", ok)
j = json.loads(subprocess.run([PERSEUS, "pcgroup", os.path.join(PCF, "g27_3.txt"), "--maxdim", "5", "--ring", "3", "--print", "json"], capture_output=True, text=True, timeout=600).stdout)
m = perseus.pc_group(**perseus.read_pcgroup(os.path.join(PCF, "g27_3.txt")), maxdim=5, ring=3)
check("command line pcgroup agrees with the module", [h["torsion"] for h in j["homology"]] == [h["torsion"] for h in m["homology"]]
      and j["cohomology"]["dims"] == m["cohomology"]["dims"] and j["cohomology"]["generators_by_degree"] == m["cohomology"]["generators_by_degree"]
      and {tuple(x[:4]): x[4] for x in j["cohomology"]["products"]} == m["cohomology"]["products"], j)
ja = json.loads(subprocess.run([PERSEUS, "pcgroup", os.path.join(PCF, "g27_3.txt"), "--maxdim", "5", "--ring", "3", "--all-products", "--print", "json"], capture_output=True, text=True, timeout=600).stdout)
ma = perseus.pc_group(**perseus.read_pcgroup(os.path.join(PCF, "g27_3.txt")), maxdim=5, ring=3, all_products=True)
check("pc_group: the full product table on request, from the command line and the module alike",
      m["cohomology"]["products_of"] == "generators" and ma["cohomology"]["products_of"] == "all basis classes" == ja["cohomology"]["products_of"]
      and {tuple(x[:4]): x[4] for x in ja["cohomology"]["products"]} == ma["cohomology"]["products"]
      and all(ma["cohomology"]["products"][k] == v for k, v in m["cohomology"]["products"].items() if k in ma["cohomology"]["products"]), (m["cohomology"]["products_of"], ma["cohomology"]["products_of"]))

SYF = os.path.join(FIX, "sylow")
r = perseus.finite_group(**perseus.read_sylow(os.path.join(SYF, "m11.txt")), maxdim=9, ring=2)
check("finite_group M11 (not solvable): H_n through Sylow subgroups", [h["torsion"] for h in r["homology"]][1:6] == [[], [], [8], [], [2]] and r["complete"], r["homology"])
c = r["cohomology"]
check("finite_group M11: mod-2 ring as in Benson-Carlson (generators in degrees 3, 4, 5)", c["dims"] == [1, 0, 0, 1, 1, 1, 1, 1, 2, 2] and c["generators_by_degree"][1:6] == [0, 0, 1, 1, 1]
      and all(len(v) == c["sylow_dims"][n] for n, b in enumerate(c["basis"]) for v in b), c)
j = json.loads(subprocess.run([PERSEUS, "sylow", os.path.join(SYF, "a5.txt"), "--maxdim", "6", "--ring", "2", "--print", "json"], capture_output=True, text=True, timeout=600).stdout)
m = perseus.finite_group(**perseus.read_sylow(os.path.join(SYF, "a5.txt")), maxdim=6, ring=2)
check("command line sylow agrees with the module", [h["torsion"] for h in j["homology"]] == [h["torsion"] for h in m["homology"]]
      and j["cohomology"]["dims"] == m["cohomology"]["dims"] and {tuple(x[:4]): x[4] for x in j["cohomology"]["products"]} == m["cohomology"]["products"], j)
ma = perseus.finite_group(**perseus.read_sylow(os.path.join(SYF, "a5.txt")), maxdim=6, ring=2, all_products=True)
check("finite_group: the full product table on request agrees with the generator products", ma["cohomology"]["products_of"] == "all basis classes"
      and m["cohomology"]["products_of"] == "generators" and ma["cohomology"]["generators_by_degree"] == m["cohomology"]["generators_by_degree"]
      and all(ma["cohomology"]["products"][k] == v for k, v in m["cohomology"]["products"].items() if k in ma["cohomology"]["products"]), ma["cohomology"]["products_of"])

PSF = os.path.join(FIX, "pcseries")
r = perseus.pc_series(**perseus.read_pcseries(os.path.join(PSF, "d64_L.txt")), maxdim=2, p=2)
check("pc_series: Ellis-King's lower central persistence matrix of D64 in degree 2",
      r["degrees"][2]["rank"] == [[3, 2, 2, 2, 2], [0, 3, 2, 2, 2], [0, 0, 3, 2, 2], [0, 0, 0, 3, 2], [0, 0, 0, 0, 3]] and (0, 4, 2) in r["degrees"][2]["barcode"], r["degrees"][2])
j = json.loads(subprocess.run([PERSEUS, "pcseries", os.path.join(PSF, "q_81_7_L.txt"), "--maxdim", "3", "--print", "json"], capture_output=True, text=True, timeout=600).stdout)
m = perseus.pc_series(**perseus.read_pcseries(os.path.join(PSF, "q_81_7_L.txt")), maxdim=3, p=3)
check("command line pcseries agrees with the module (ranks, bars, maps)", j["p"] == 3 and all(
      a["rank"] == b["rank"] and [list(x) for x in b["barcode"]] == a["barcode"] and a["maps"] == b["maps"] for a, b in zip(j["degrees"], m["degrees"])), j["degrees"][1])

# ---- torsion trajectories and classification by Howell sweeps (torsion.hpp) ---------------------------------
def smith_torsion(cols, nrows):
    """invariant factors > 1 of the integer matrix with these columns (lists of (row, coef)), by elimination"""
    A = [[0] * len(cols) for _ in range(nrows)]
    for j, c in enumerate(cols):
        for r, v in c: A[r][j] += v
    out = []
    while A and A[0]:
        nz = [(abs(A[i][j]), i, j) for i in range(len(A)) for j in range(len(A[0])) if A[i][j]]
        if not nz: break
        _, i, j = min(nz)
        A[0], A[i] = A[i], A[0]
        for row in A: row[0], row[j] = row[j], row[0]
        done = False
        while not done:
            done = True
            p = A[0][0]
            for i in range(1, len(A)):
                q = A[i][0] // p
                if q: A[i] = [a - q * b for a, b in zip(A[i], A[0])]
                if A[i][0]: done = False
            for j in range(1, len(A[0])):
                q = A[0][j] // p
                if q:
                    for row in A: row[j] -= q * row[0]
                if A[0][j]: done = False
            if not done:
                nz = [(abs(A[i][0]), i, 0) for i in range(len(A)) if A[i][0]] + [(abs(A[0][j]), 0, j) for j in range(len(A[0])) if A[0][j]]
                _, i, j = min(nz)
                A[0], A[i] = A[i], A[0]
                for row in A: row[0], row[j] = row[j], row[0]
            elif any(A[i][j] % p for i in range(1, len(A)) for j in range(1, len(A[0]))):
                i = next(i for i in range(1, len(A)) if any(x % p for x in A[i][1:]))
                A[0] = [a + b for a, b in zip(A[0], A[i])]; done = False
        out.append(abs(A[0][0]))
        A = [row[1:] for row in A[1:]]
    return sorted(x for x in out if x > 1)

# a hand example: e1, e2 at 0; 4 e1 at 1; 6 e2 at 2; 2 e1 + 3 e2 at 3; e1 at 4
_hand = [(1, 0.0, []), (1, 0.0, []), (2, 1.0, [(0, 4)]), (2, 2.0, [(1, 6)]), (2, 3.0, [(0, 2), (1, 3)]), (2, 4.0, [(0, 1)])]
_t = perseus.cells(_hand, torsion=True)["torsion"]
check("torsion: hand example Z/4, Z/2 + Z/12, Z/12, Z/3",
      [(s["value"], s["invariant_factors"]) for s in _t["degrees"][1]] == [(1.0, [4]), (2.0, [2, 12]), (3.0, [12]), (4.0, [3])] and _t["proof"], _t)
_t2 = perseus.cells(_hand, torsion=True, reduce=False)["torsion"]
check("torsion: the same without the algebraic Morse reduction", _t2["degrees"] == _t["degrees"])
# random filtered chain complexes (n 1-cells, m 2-cells), tied and distinct values, against a Smith form at every t
_rng = random.Random(7); _bad = 0; _seen = 0
for _it in range(12):
    _n, _m, _tied = 8, 16, _it % 2 == 0
    _b1 = [0.0 if _tied else round(_rng.random(), 6) for _ in range(_n)]
    _cells = [(1, b, []) for b in _b1]; _cols = []
    for _j in range(_m):
        _bd = [(f, _rng.choice([1, -1, 2, -2, 3, 4, 6, 8, 9, 12])) for f in _rng.sample(range(_n), _rng.choice([1, 2, 2, 3]))]
        _v = float(_j + 1) if _tied else round(max(_b1[f] for f, _ in _bd) + _rng.random() + 1e-6, 6)
        _cells.append((2, _v, _bd)); _cols.append((_v, _bd))
    _r = perseus.cells(_cells, torsion=True)
    _st = _r["torsion"]["degrees"].get(1, [])
    for _v in sorted({v for v, _ in _cols}):
        _want = smith_torsion([bd for v, bd in _cols if v <= _v], _n)
        _got = []
        for s in _st:
            if s["value"] <= _v: _got = s["invariant_factors"]
        _seen += bool(_want); _bad += _got != _want
check("torsion: random chain complexes, every t against a Smith form", _bad == 0 and _seen > 20, (_bad, _seen))
# Linial-Meshulam n = 60 (above): sweeps and reruns agree, and the torsion's primes are exactly the bad primes
_rs = perseus.simplices(_lm, classify_by="sweep", torsion=True); _rr = perseus.simplices(_lm, classify_by="rerun")
_same = sorted(_rs["certificate"]["bad_primes"]) == sorted(_rr["certificate"]["bad_primes"]) and all(
    all(len(_rs["certificate"]["bad_primes"][p][k]) == len(_rr["certificate"]["bad_primes"][p][k]) and
        all(np.array_equal(x, y) for x, y in zip(_rs["certificate"]["bad_primes"][p][k], _rr["certificate"]["bad_primes"][p][k]))
        for k in ("only_over_Fp", "only_over_Q")) for p in _rr["certificate"]["bad_primes"])
check("classification: sweeps and reruns give the same bad primes and intervals (LM n = 60)", _same)
_tp = sorted({p for s in _rs["torsion"]["degrees"][1] for p in s["primary"]})
check("torsion: on a Linial-Meshulam process its primes are the bad primes (Theorem 4)", _tp == sorted(_rr["certificate"]["bad_primes"]), (_tp, sorted(_rr["certificate"]["bad_primes"])))
_p = subprocess.run([PERSEUS, "simplices", _path, "--print", "json", "--torsion"], capture_output=True, text=True, timeout=600); _j = json.loads(_p.stdout)
check("torsion: command line agrees with the module",
      [(s["value"], [int(x) for x in s["invariant_factors"]]) for s in _j["torsion"]["degrees"][0]["steps"]] ==
      [(s["value"], s["invariant_factors"]) for s in _rs["torsion"]["degrees"][1]] and _j["torsion"]["proof"] is True)
_rb = perseus.cells(read_cells(os.path.join(FIX, "crt_alias_1.cells")), torsion=True)["torsion"]["degrees"][1]
check("torsion: 98-bit prime with its exponent", _rb[-1]["primary"].get(294719631262177595678027596519) == [1] and
      _rb[-1]["invariant_factors"] == [2**2 * 3 * 5**2 * 7**2 * 4909 * 294719631262177595678027596519], _rb)
_rg = perseus.cells(read_cells(os.path.join(FIX, "crt_group.cells")), torsion=True)["torsion"]
check("torsion: unfactored composite reported as a group", _rg["groups"] == [36028797018976327 * 1152921504606848017] and _rg["complete"], _rg)

print(f"\n{results['pass']} passed, {results['fail']} failed")
for f in failures[:30]: print("  ", f)
sys.exit(1 if results["fail"] else 0)
