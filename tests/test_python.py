#!/usr/bin/env python3
"""The Python module against the command line on the same inputs: barcodes, bad primes and certificate
fields, in every certificate mode. Needs the module installed (cd python && pip install .) and ./perseus."""
import json, os, subprocess, sys, tempfile, random
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

print(f"\n{results['pass']} passed, {results['fail']} failed")
for f in failures[:30]: print("  ", f)
sys.exit(1 if results["fail"] else 0)
