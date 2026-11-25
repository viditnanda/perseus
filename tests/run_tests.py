#!/usr/bin/env python3
"""Differential tests for Perseus 6.

  A. Rips engine vs Ripser (exact interval multisets) over F_p, and the Q barcode vs Ripser mod a good prime
  B. Rips engine vs the explicit engine on the same complexes: barcodes and certificates (bad primes)
  C. explicit engine vs the independent Python reference p2ref (Z-elimination + field reductions)
  D. old Perseus 4/5 files, converted by tools/convert_legacy.py, vs the Perseus 5 binary (F_2)
  E. flag / cubical input routes against each other
  F. certificate regressions: inputs whose pivots vanish or look like units modulo fixed primes, and
     unfactorable pivot values (exact expected answers, in auto, fast and proved mode)
  G. implicit cubical engine vs the explicit complex; sparse Rips route vs the implicit Rips engine

Paths of external oracles are taken from environment variables (RIPSER, PERSEUS5, P2REF); a section
whose oracle is missing is skipped.
"""
import itertools, json, os, random, subprocess, sys, tempfile, math
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
PERSEUS = os.environ.get("PERSEUS", os.path.join(HERE, "..", "perseus"))
RIPSER = os.environ.get("RIPSER", "/home/claude/ext/ripser/ripser-coeff")
PERSEUS5 = os.environ.get("PERSEUS5", "/home/claude/build/perseus")
P2REF = os.environ.get("P2REF", "/home/claude/perseus2")
BENCH = os.environ.get("BENCH", "/home/claude/bench/bin")
TMP = tempfile.mkdtemp(prefix="p6test_")
results = Counter(); failures = []

def check(name, ok, detail=""):
    results["pass" if ok else "fail"] += 1
    if not ok:
        failures.append(name + (": " + detail if detail else ""))
        print("FAIL", name, detail[:300], flush=True)

def run(cmd, timeout=600):
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(f"{cmd} failed: {r.stderr[:500]}")
    return r.stdout, r.stderr

def key(b, d):
    return ("%.6g" % b, "inf" if d is None or d == float("inf") else "%.6g" % d)

def parse_ripser(out):
    bars = Counter(); dim = None
    for line in out.splitlines():
        if line.startswith("persistence intervals in dim"):
            dim = int(line.split()[-1].rstrip(":")); continue
        line = line.strip()
        if dim is None or not line.startswith("["): continue
        a, b = line[1:].rstrip(")").split(",")
        bars[(dim,) + key(float(a), float(b) if b.strip() else None)] += 1
    return bars

def perseus_json(args):
    out, err = run([PERSEUS] + args + ["--print", "json"])
    j = json.loads(out)
    bars = Counter()
    for k, lst in enumerate(j["barcode"]):
        for b, d in lst:
            bars[(k,) + key(b, d)] += 1
    return bars, j

def perseus_field(args, p):
    out, _ = run([PERSEUS] + args + ["--coeff", str(p), "--print", "ripser"])
    return parse_ripser(out)

def restrict(bars, maxdim):
    return Counter({k: v for k, v in bars.items() if k[0] <= maxdim})

# ---------------------------------------------------------------- A. vs Ripser
def section_A(only=None):
    if not os.path.exists(RIPSER): print("skip A (no ripser)"); return
    cases = [("sphere", 2, None), ("o3_1024", 3, 1.8), ("dragon", 1, None), ("fractal", 2, None),
             ("random16", 4, None), ("randdist120", 2, None), ("randdist180", 1, None), ("randdist250", 1, None)]
    for f, dim, thr in cases:
        if only and f not in only: continue
        path = f"{BENCH}/{f}.bin"
        if not os.path.exists(path): continue
        ta = ["--threshold", str(thr)] if thr else []
        q, j = perseus_json(["bin32", path, "--maxdim", str(dim)] + ta)
        bad = [b["p"] for b in j["bad_primes"]]
        for p in (2, 3, 5, 7, 127):
            rp = parse_ripser(run([RIPSER, "--format", "binary", "--dim", str(dim), "--modulus", str(p)] + ta + [path])[0])
            ours = perseus_field(["bin32", path, "--maxdim", str(dim)] + ta, p)
            check(f"A {f} H<={dim} F_{p} vs ripser", ours == rp, f"diff ours-rip={dict(ours-rp)} rip-ours={dict(rp-ours)}")
            if p not in bad and p != 127 and p in j["candidates"]:
                pass
            if p not in bad:
                check(f"A {f} H<={dim} Q vs ripser F_{p} (p good)", q == rp, f"{dict(q-rp)} {dict(rp-q)}")
            else:
                check(f"A {f} H<={dim} Q differs from F_{p} (p bad)", q != rp)
        print(f"A {f}: bad primes {bad}", flush=True)

# ---------------------------------------------------------------- B. Rips engine vs explicit engine
def write_lower(D, path):
    n = len(D)
    with open(path, "w") as fh:
        for i in range(1, n):
            fh.write(",".join("%.17g" % D[i][j] for j in range(i)) + "\n")

def section_B(trials=60):
    rnd = random.Random(1)
    for t in range(trials):
        n = rnd.randint(6, 40); dim = rnd.choice([1, 2, 2, 3])
        kind = rnd.choice(["uniform", "integer", "euclid"])
        if kind == "euclid":
            dd_ = rnd.choice([2, 3, 4]); P = [[rnd.random() for _ in range(dd_)] for _ in range(n)]
            D = [[math.dist(P[i], P[j]) for j in range(n)] for i in range(n)]
        else:
            D = [[0.0] * n for _ in range(n)]
            for i in range(n):
                for j in range(i):
                    D[i][j] = D[j][i] = (rnd.random() if kind == "uniform" else float(rnd.randint(1, 6)))
        path = f"{TMP}/b{t}.txt"; write_lower(D, path)
        thr = ["--threshold", "1e9"]
        a, ja = perseus_json(["lower", path, "--maxdim", str(dim)] + thr)
        b, jb = perseus_json(["lower", path, "--maxdim", str(dim), "--engine", "explicit"] + thr)
        check(f"B{t} n={n} {kind} H<={dim} Q barcode rips==explicit", a == b, f"{dict(a-b)} {dict(b-a)}")
        ba = sorted(x["p"] for x in ja["bad_primes"]); bb = sorted(x["p"] for x in jb["bad_primes"])
        check(f"B{t} bad primes rips==explicit", ba == bb, f"{ba} vs {bb}")
        cc, jc = perseus_json(["lower", path, "--maxdim", str(dim), "--collapse", "full"] + thr)
        check(f"B{t} Q barcode collapsed rips==explicit", cc == b, f"{dict(cc-b)} {dict(b-cc)}")
        bc = sorted(x["p"] for x in jc["bad_primes"])
        check(f"B{t} bad primes collapsed==explicit", bc == bb, f"{bc} vs {bb}")
        for p in (2, 3):
            fa = perseus_field(["lower", path, "--maxdim", str(dim)] + thr, p)
            fb = perseus_field(["lower", path, "--maxdim", str(dim), "--engine", "explicit"] + thr, p)
            check(f"B{t} F_{p} rips==explicit", fa == fb, f"{dict(fa-fb)} {dict(fb-fa)}")

# ---------------------------------------------------------------- C. explicit engine vs p2ref
def to_cells(C, path):
    with open(path, "w") as fh:
        for c in range(len(C)):
            fh.write(f"{C.dim[c]} {C.birth[c]!r}" + "".join(f" {s} {v}" for s, v in sorted(C.bd[c].items())) + "\n")

def section_C(trials=80):
    if not os.path.exists(os.path.join(P2REF, "p2ref.py")): print("skip C (no p2ref)"); return
    sys.path.insert(0, P2REF); import p2ref
    rnd = random.Random(2)
    for t in range(trials):
        kind = rnd.choice(["lm", "simp", "chain"])
        if kind == "lm":     # Linial-Meshulam style: all edges at 0, random triangles in random order
            n = rnd.randint(5, 9)
            from itertools import product, combinations
            tris = list(combinations(range(n), 3)); rnd.shuffle(tris); tris = tris[:rnd.randint(1, len(tris))]
            order = {tr: i + 1 for i, tr in enumerate(tris)}
            facets = [e for e in combinations(range(n), 2)] + tris
            C, _ = p2ref.simplicial(facets, birth_of=lambda s: order.get(s, 0) if len(s) == 3 else 0)
        elif kind == "simp":
            n = rnd.randint(5, 8)
            facets = [tuple(sorted(rnd.sample(range(n), rnd.randint(2, 4)))) for _ in range(rnd.randint(3, 12))]
            bt = {}
            C, _ = p2ref.simplicial(facets, birth_of=lambda s: bt.setdefault(s, rnd.randint(0, 5)))
        else:                # torsion-rich chain complex: one vertex, cycles, random 2-cells
            C = p2ref.FComplex(); v = C.add(0, 0)
            m1 = rnd.randint(1, 4); e = [C.add(1, rnd.randint(0, 3)) for _ in range(m1)]
            for _ in range(rnd.randint(1, 5)):
                bd = {x: rnd.choice([0, 1, -1, 2, -2, 3, 4, 6, -5, 7, 9, 12]) for x in e}
                C.add(2, max(C.birth[x] for x in e) + rnd.randint(0, 3), bd)
        C.check()
        path = f"{TMP}/c{t}.txt"; to_cells(C, path)
        ref = p2ref.analyse(C)
        ours, j = perseus_json(["cells", path])
        refq = Counter((d,) + key(b, e) for d, b, e in ref["Q"])
        check(f"C{t} {kind} Q barcode vs p2ref", ours == refq, f"{dict(ours-refq)} {dict(refq-ours)}")
        ob = sorted(x["p"] for x in j["bad_primes"])
        check(f"C{t} {kind} bad primes vs p2ref", ob == sorted(ref["bad"]), f"{ob} vs {ref['bad']}")
        for p in ref["Fp"]:
            refp = Counter((d,) + key(b, e) for d, b, e in ref["Fp"][p])
            op = perseus_field(["cells", path], p)
            check(f"C{t} {kind} F_{p} vs p2ref", op == refp, f"{dict(op-refp)} {dict(refp-op)}")

# ---------------------------------------------------------------- D. legacy formats vs Perseus 5
def read_legacy_out(prefix, maxdim=12):
    bars = Counter()
    for k in range(maxdim + 1):
        fn = f"{prefix}_{k}.txt"
        if not os.path.exists(fn): continue
        for line in open(fn):
            t = line.split()
            if len(t) == 2:
                d = float(t[1]); bars[(k,) + key(float(t[0]), None if d == -1 else d)] += 1
    return bars

CONVERT = os.path.join(HERE, "..", "tools", "convert_legacy.py")
def legacy_pair(fmt, text, name):
    """old file -> tools/convert_legacy.py -> Perseus 6, against the Perseus 5 executable on the old file"""
    path = f"{TMP}/{name}.txt"
    open(path, "w").write(text)
    p5 = f"{TMP}/{name}_p5"; p6 = f"{TMP}/{name}_p6"
    subprocess.run([PERSEUS5, fmt, path, p5], capture_output=True, timeout=120, stdin=subprocess.DEVNULL)
    sys.path.insert(0, os.path.dirname(CONVERT)); import convert_legacy
    args = convert_legacy.convert(fmt, path, f"{TMP}/{name}_new.txt")
    run([PERSEUS] + args[:2] + [p6] + args[2:])
    maxdim = int(args[args.index("--maxdim") + 1]) if "--maxdim" in args else 12
    a, b = restrict(read_legacy_out(p5), maxdim), read_legacy_out(p6)
    check(f"D {name} {fmt} (converted) vs Perseus5", a == b, f"p5-p6={dict(a-b)} p6-p5={dict(b-a)}")

def section_D(trials=25):
    if not os.path.exists(PERSEUS5): print("skip D (no perseus5)"); return
    rnd = random.Random(3)
    for t in range(trials):
        # nmfsimtop: ptdim 1, simplices of varying dimension with integer births >= 1
        lines = ["1"]
        for _ in range(rnd.randint(3, 15)):
            d = rnd.randint(1, 3); vs = rnd.sample(range(1, 9), d + 1)
            lines.append(" ".join(map(str, [d] + vs + [rnd.randint(1, 6)])))
        legacy_pair("nmfsimtop", "\n".join(lines) + "\n", f"nmf{t}")
        # simtop: uniform dimension 2
        lines = ["2", "1"]
        for _ in range(rnd.randint(3, 15)):
            vs = rnd.sample(range(1, 9), 3); lines.append(" ".join(map(str, vs + [rnd.randint(1, 6)])))
        legacy_pair("simtop", "\n".join(lines) + "\n", f"sim{t}")
        # cubtop: random 2d or 3d grid, some periodic axes, births in 1..5 or -1 (absent)
        d = rnd.choice([2, 2, 3])
        ext = [rnd.randint(1, 5 if d == 2 else 3) for _ in range(d)]
        sext = [(-e if rnd.random() < 0.3 else e) for e in ext]
        tot = 1
        for e in ext: tot *= e
        bs = [rnd.choice([-1, 1, 2, 3, 4, 5]) for _ in range(tot)]
        legacy_pair("cubtop", "\n".join(map(str, [d] + sext + bs)) + "\n", f"cub{t}")
        # scubtop
        lines = [str(d)]
        for _ in range(rnd.randint(2, 12)):
            lines.append(" ".join(map(str, [rnd.randint(0, 3) for _ in range(d)] + [rnd.randint(1, 5)])))
        legacy_pair("scubtop", "\n".join(lines) + "\n", f"scub{t}")
        # distmat
        n = rnd.randint(4, 9)
        M = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i):
                M[i][j] = M[j][i] = round(rnd.random() * 10, 3)
        g, s, N, dd = 0.5, round(rnd.choice([0.5, 1.0, 1.5]), 2), rnd.randint(3, 12), rnd.choice([2, 3])
        txt = f"{n}\n{g} {s} {N} {dd}\n" + "\n".join(" ".join(map(str, r)) for r in M) + "\n"
        legacy_pair("distmat", txt, f"dm{t}")
        # rips (ball growing): 2d points with radii
        pts = [(round(rnd.random() * 4, 3), round(rnd.random() * 4, 3), round(rnd.random() * 0.5, 2)) for _ in range(rnd.randint(4, 9))]
        txt = f"2\n{rnd.choice([0.5, 1.0])} {rnd.choice([0.1, 0.25])} {rnd.randint(4, 15)}\n" + "\n".join(" ".join(map(str, q)) for q in pts) + "\n"
        legacy_pair("brips", txt, f"br{t}")
        # rips (vertex births, fixed radii)
        pts = [(round(rnd.random() * 4, 3), round(rnd.random() * 4, 3), round(0.3 + rnd.random(), 2), rnd.randint(1, 6)) for _ in range(rnd.randint(4, 9))]
        legacy_pair("rips", "2\n" + "\n".join(" ".join(map(str, q)) for q in pts) + "\n", f"vr{t}")
        # cellcomp: a simplicial complex written as cells (coefficients +-1), random births
        from itertools import combinations
        n = rnd.randint(4, 7)
        facets = [tuple(sorted(rnd.sample(range(n), rnd.randint(2, 3)))) for _ in range(rnd.randint(2, 8))]
        simp = sorted({s for f in facets for k in range(1, len(f) + 1) for s in combinations(f, k)}, key=lambda s: (len(s), s))
        bt = {}
        for s in simp:
            b = rnd.randint(1, 4)
            for i in range(len(s) if len(s) > 1 else 0):
                b = max(b, bt[s[:i] + s[i + 1:]])
            bt[s] = b
        lines = []; idx = {}
        for k in range(1, 4):
            layer = [s for s in simp if len(s) == k]
            if not layer: break
            lines.append("-1")
            for j, s in enumerate(layer):
                idx[s] = j
                if k == 1: lines.append(f"{bt[s]} 0")
                else:
                    fs = [s[:i] + s[i + 1:] for i in range(k)]
                    lines.append(f"{bt[s]} {k} " + " ".join(f"{(-1) ** i} {idx[f]}" for i, f in enumerate(fs)))
        # Perseus 5's own cellcomp reader is broken (it mangles a filled triangle), so the converted
        # cellcomp file is checked against the same complex written directly in the 'cells' format
        path = f"{TMP}/cc{t}.txt"; open(path, "w").write("\n".join(lines) + "\n")
        glob = {}; cl = []
        for s in simp:
            glob[s] = len(cl)
            fs = [s[:i] + s[i + 1:] for i in range(len(s))] if len(s) > 1 else []
            cl.append(f"{len(s) - 1} {bt[s]}" + "".join(f" {glob[f]} {(-1) ** i}" for i, f in enumerate(fs)))
        path2 = f"{TMP}/cc{t}_cells.txt"; open(path2, "w").write("\n".join(cl) + "\n")
        sys.path.insert(0, os.path.dirname(CONVERT)); import convert_legacy
        args = convert_legacy.convert("cellcomp", path, f"{TMP}/cc{t}_conv.txt")
        a = perseus_field(args[:2], 2); b = perseus_field(["cells", path2], 2)
        check(f"D cc{t} cellcomp (converted) == cells", a == b, f"{dict(a-b)} {dict(b-a)}")


# ---------------------------------------------------------------- E. new formats vs independent constructions
def section_E(trials=40):
    rnd = random.Random(5)
    from itertools import combinations, product
    for t in range(trials):
        # flag with vertex births vs explicitly enumerated cliques
        n = rnd.randint(5, 14); vb = [rnd.choice([0, 0, 1, 2, 3]) for _ in range(n)]
        E = {}
        for i, j in combinations(range(n), 2):
            if rnd.random() < 0.6: E[(i, j)] = max(vb[i], vb[j]) + rnd.randint(0, 6)
        fpath = f"{TMP}/fl{t}.txt"
        with open(fpath, "w") as fh:
            for i in range(n): fh.write(f"{i} {vb[i]}\n")
            for (i, j), w in E.items(): fh.write(f"{i} {j} {w}\n")
        adj = {i: set() for i in range(n)}
        for (i, j) in E: adj[i].add(j); adj[j].add(i)
        spath = f"{TMP}/fs{t}.txt"; md = rnd.choice([1, 2])
        with open(spath, "w") as fh:
            for k in range(1, md + 3):
                for s in combinations(range(n), k):
                    if all(b in adj[a] for a, b in combinations(s, 2)):
                        bt = max([vb[v] for v in s] + [E[(a, b)] for a, b in combinations(s, 2)])
                        fh.write(" ".join(map(str, s)) + f" {bt}\n")
        a, ja = perseus_json(["flag", fpath, "--maxdim", str(md)])
        b, jb = perseus_json(["simplices", spath, "--maxdim", str(md)])
        check(f"E{t} flag(vertex births) == simplices", a == b, f"{dict(a-b)} {dict(b-a)}")
        check(f"E{t} flag bad primes == simplices", sorted(x["p"] for x in ja["bad_primes"]) == sorted(x["p"] for x in jb["bad_primes"]))
        # flag without vertex births vs lower-distance with inf for missing edges
        fpath2 = f"{TMP}/fl0_{t}.txt"; lpath = f"{TMP}/ld{t}.txt"
        with open(fpath2, "w") as fh:
            for i in range(n): fh.write(f"{i} 0\n")
            for (i, j), w in E.items(): fh.write(f"{i} {j} {w}\n")
        with open(lpath, "w") as fh:
            for i in range(1, n): fh.write(",".join(str(E.get((j, i), "inf")) for j in range(i)) + "\n")
        a = perseus_field(["flag", fpath2, "--maxdim", str(md)], 2)
        b = perseus_field(["lower-distance", lpath, "--maxdim", str(md), "--threshold", "1e9"], 2)
        check(f"E{t} flag == lower-distance with inf", a == b, f"{dict(a-b)} {dict(b-a)}")
        # sparse cubes vs dense cubical with absent cells
        d = rnd.choice([2, 3]); shape = [rnd.randint(1, 5) for _ in range(d)]
        import numpy as np
        A = np.where(np.random.default_rng(t).random(shape) < 0.25, np.nan, np.random.default_rng(t + 99).integers(1, 9, size=shape).astype(float))
        cpath = f"{TMP}/cu{t}.txt"; spath2 = f"{TMP}/cs{t}.txt"
        with open(cpath, "w") as fh: fh.write(" ".join(map(str, shape)) + "\n" + "\n".join("nan" if np.isnan(x) else str(int(x)) for x in A.flatten()) + "\n")
        with open(spath2, "w") as fh:
            for idx in np.ndindex(*shape):
                if not np.isnan(A[idx]): fh.write(" ".join(map(str, idx)) + f" {int(A[idx])}\n")
        if np.isnan(A).all(): continue
        a = perseus_field(["cubical", cpath], 3); b = perseus_field(["cubes", spath2], 3)
        check(f"E{t} cubical == cubes", a == b, f"{dict(a-b)} {dict(b-a)}")
        # cubes --integral --image (algebraic Morse reduction, packed builder) vs the same complex
        # written out cell by cell by an independent builder and read as cells
        births = {}
        for idx in np.ndindex(*shape):
            if np.isnan(A[idx]): continue
            for off in product(range(3), repeat=d):
                c = tuple(2 * i + o for i, o in zip(idx, off))
                births[c] = min(births.get(c, 99), int(A[idx]))
        order = sorted(births, key=lambda c: (sum(x & 1 for x in c), c))
        ids = {c: k for k, c in enumerate(order)}
        cellpath = f"{TMP}/cc{t}.txt"
        with open(cellpath, "w") as fh:
            for c in order:
                row, j = [], 0
                for ax in range(d):
                    if c[ax] & 1:
                        s = -1 if j & 1 else 1
                        for dr in (1, -1):
                            f = list(c); f[ax] += dr; row += [ids[tuple(f)], s * dr]
                        j += 1
                fh.write(" ".join(map(str, [sum(x & 1 for x in c), births[c]] + row)) + "\n")
        lo_t, hi_t = sorted(rnd.sample(range(1, 10), 2))
        _, ja = perseus_json(["cubes", spath2, "--integral", "--image", f"{lo_t},{hi_t}"])
        _, jb = perseus_json(["cells", cellpath, "--integral", "--image", f"{lo_t},{hi_t}"])
        check(f"E{t} cubes integral image == cells", ja.get("integral_image") == jb.get("integral_image") and "integral_image" in ja,
              f"{ja.get('integral_image')} {jb.get('integral_image')} {ja.get('integral_error')}")


# ---------------------------------------------------------------- F. certificate regressions
FIX = os.path.join(HERE, "fixtures")
M4 = [2147483579, 2147483587, 2147483629, 2147483647]
M1 = [2, 3, 5, 7, 4909, 294719631262177595678027596519]
def section_F():
    cases = [("cells", "crt_alias_0.cells", M4), ("cells", "crt_alias_1.cells", M1),
             ("cells", "crt_unit_0.cells", M4), ("cells", "crt_unit_1.cells", M1),
             ("simplices", "crt_unit_0.simplices", M4), ("simplices", "crt_unit_1.simplices", M1),
             ("cells", "crt_group_split.cells", [36028797018976327, 144115188075855881, 288230376151711849, 1152921504606848017])]
    for fmt, f, expect in cases:
        for mode in ("auto", "fast", "proved"):
            bars, j = perseus_json([fmt, os.path.join(FIX, f), "--mode", mode])
            got = sorted(int(b["p"]) for b in j["bad_primes"])
            check(f"F {f} {mode} bad primes", got == sorted(expect), f"{got}")
            check(f"F {f} {mode} complete", j["complete"] and j["guarantee"]["barcode_certified"])
            # the rational barcode: H1 bars all finite, nothing in H2 (the determinant is nonzero)
            check(f"F {f} {mode} no essential H1/H2", not any(k >= 1 and d == "inf" for (k, b, d) in bars), f"{[x for x in bars if x[0] >= 1 and x[2] == 'inf']}")
            check(f"F {f} {mode} H2 empty", not any(k == 2 for (k, b, d) in bars))
    for mode in ("auto", "fast", "proved"):
        _, j = perseus_json(["cells", os.path.join(FIX, "crt_group.cells"), "--mode", mode])
        g = [b for b in j["bad_primes"] if b.get("every_prime_factor_of")]
        check(f"F crt_group {mode}: one group entry for p*q", len(g) == 1 and int(g[0]["p"]) == 36028797018976327 * 1152921504606848017, f"{j['bad_primes']}")
    # a composite cofactor N appears in degree 1 (inside 3N) and again in degree 2, where it changes no bar
    # of positive length: N must keep its earliest degree, or its primes are missed
    N1 = 2305843009213706309 * 2305843010201348287; N2 = 2306033382932562277 * 2306753820319319929
    for f, expect in (("crt_group_degree.cells", [(3, False), (N1, True)]), ("crt_shared_degree.cells", [(2, False), (3, False), (N2, True)])):
        for mode in ("auto", "fast", "proved"):
            _, j = perseus_json(["cells", os.path.join(FIX, f), "--mode", mode])
            got = sorted((int(b["p"]), bool(b.get("every_prime_factor_of"))) for b in j["bad_primes"])
            check(f"F {f} {mode} bad primes", got == sorted(expect) and j["complete"], f"{got}")
            check(f"F {f} {mode} groups from degree 1", all(b["first_degree"] == 1 for b in j["bad_primes"]), f"{j['bad_primes']}")
    # fast mode: a fresh run draws a 256-bit seed and certifies; a supplied seed reproduces the run exactly
    # but certifies nothing (the bound is a probability over a seed drawn after the input is fixed)
    a = perseus_json(["cells", os.path.join(FIX, "crt_alias_1.cells"), "--mode", "fast", "--seed", "1"])[1]
    a2 = perseus_json(["cells", os.path.join(FIX, "crt_alias_1.cells"), "--mode", "fast", "--seed", "0x1"])[1]
    b = perseus_json(["cells", os.path.join(FIX, "crt_alias_1.cells"), "--mode", "fast", "--seed", "2"])[1]
    f1 = perseus_json(["cells", os.path.join(FIX, "crt_alias_1.cells"), "--mode", "fast"])[1]
    check("F seeds: same answer", a["barcode"] == b["barcode"] == f1["barcode"] and a["bad_primes"] == b["bad_primes"] == f1["bad_primes"])
    check("F seeds: different primes", a["guarantee"]["moduli"] != b["guarantee"]["moduli"])
    check("F seeds: a seed reproduces its run", a["guarantee"]["moduli"] == a2["guarantee"]["moduli"])
    check("F seeds: supplied seed is not certified", not a["guarantee"]["barcode_certified"] and not a["complete"] and a["guarantee"]["seed_supplied"]
          and a["bad_primes_status"].startswith("conditional"), f"{a['guarantee']} {a['complete']}")
    check("F seeds: fresh seed is 256 bits and certifies", len(f1["guarantee"]["seed"]) == 64 and f1["guarantee"]["barcode_certified"] and f1["complete"]
          and not f1["guarantee"]["seed_supplied"])
    r = perseus_json(["cells", os.path.join(FIX, "crt_alias_1.cells"), "--mode", "fast", "--seed", f1["guarantee"]["seed"]])[1]
    check("F seeds: the reported seed reproduces a fresh run", r["guarantee"]["moduli"] == f1["guarantee"]["moduli"])

# ---------------------------------------------------------------- G. new engines vs explicit references
def section_G(trials=60):
    import numpy as np
    rnd = random.Random(5); rng = np.random.default_rng(5)
    for t in range(trials):
        d = rnd.randint(1, 4); shape = tuple(rnd.randint(1, 6 if d < 4 else 3) for _ in range(d))
        A = rng.integers(0, rnd.randint(2, 10), size=shape).astype(float)
        if rnd.random() < 0.3: A[rng.random(shape) < 0.15] = np.nan
        path = f"{TMP}/g{t}.npy"; np.save(path, A)
        per = [str(i) for i in range(d) if rnd.random() < 0.35]
        pa = ["--periodic", ",".join(per)] if per else []
        for p in (2, 3):
            a = perseus_field(["npy", path] + pa, p); b = perseus_field(["npy", path, "--engine", "explicit"] + pa, p)
            check(f"G{t} cubical {shape} per={per} F_{p} implicit == explicit", a == b, f"{dict(a-b)} {dict(b-a)}")
        a, ja = perseus_json(["npy", path] + pa); b, jb = perseus_json(["npy", path, "--engine", "explicit"] + pa)
        check(f"G{t} cubical Q implicit == explicit", a == b and sorted(str(x["p"]) for x in ja["bad_primes"]) == sorted(str(x["p"]) for x in jb["bad_primes"]))
        # the classic implicit engine and the discrete-gradient engine (default where it applies), with zero-length bars
        c, jc = perseus_json(["npy", path, "--engine", "classic"] + pa)
        check(f"G{t} cubical Q classic == default", c == a and sorted(str(x["p"]) for x in jc["bad_primes"]) == sorted(str(x["p"]) for x in ja["bad_primes"]))
        z1 = perseus_field(["npy", path, "--zero-length"] + pa, 2); z2 = perseus_field(["npy", path, "--zero-length", "--engine", "explicit"] + pa, 2)
        check(f"G{t} cubical F_2 zero-length implicit == explicit", z1 == z2, f"{dict(z1-z2)} {dict(z2-z1)}")
    # the discrete-gradient engine on 4D white noise and quantised fields large enough for long V-paths and
    # Morse path sums of +-2 (against the explicit complex, over F_2, F_3 and Q)
    for t, (n, q) in enumerate([(9, None), (10, None), (10, 4), (8, 2)]):
        A = rng.random((n,) * 4)
        if q: A = np.floor(A * q)
        path = f"{TMP}/g4_{t}.npy"; np.save(path, A)
        for p in (2, 3):
            a = perseus_field(["npy", path], p); b = perseus_field(["npy", path, "--engine", "explicit"], p)
            check(f"G 4D {n}^4 levels={q} F_{p} morse == explicit", a == b, f"{dict(a-b)} {dict(b-a)}")
        a, ja = perseus_json(["npy", path]); b, jb = perseus_json(["npy", path, "--engine", "explicit"])
        check(f"G 4D {n}^4 levels={q} Q morse == explicit", a == b and sorted(str(x["p"]) for x in ja["bad_primes"]) == sorted(str(x["p"]) for x in jb["bad_primes"]))
    for t in range(12):
        n = rnd.randint(30, 300); dd = rnd.choice([2, 3, 4]); md = rnd.choice([1, 2])
        P = rng.random((n, dd)); path = f"{TMP}/sp{t}.txt"; np.savetxt(path, P)
        thr = str(rnd.choice([0.15, 0.25, 0.4]))
        for p in (2, "Q"):
            args = ["points", path, "--maxdim", str(md), "--threshold", thr]
            if p == "Q":
                a, ja = perseus_json(args + ["--engine", "sparse"]); b, jb = perseus_json(args + ["--engine", "implicit"])
            else:
                a = perseus_field(args + ["--engine", "sparse"], p); b = perseus_field(args + ["--engine", "implicit"], p)
            check(f"G sparse rips n={n} d={dd} H<={md} {p} sparse == implicit", a == b, f"{dict(a-b)} {dict(b-a)}")
    # the grid must not lose edges when coordinates differ hugely in size (cells of nearby points round apart)
    for t, (big, near) in enumerate([(-2.0**55, (3, 5)), (2.0**60, (-1, 1.5)), (-2.0**45, (0.1, 2.9)), (1e30, (7, 9))]):
        path = f"{TMP}/bigc{t}.txt"
        with open(path, "w") as fh:
            fh.write(f"{big!r} 0 0\n{near[0]} 0 0\n{near[1]} 0 0\n")
            for i in range(20): fh.write(f"{near[0] + 0.137 * i} {0.01 * i} 0\n")
        a = perseus_field(["points", path, "--threshold", "3", "--engine", "sparse"], 2)
        b = perseus_field(["points", path, "--threshold", "3", "--engine", "implicit"], 2)
        check(f"G sparse grid, large coordinates {big:g}", a == b, f"{dict(a-b)} {dict(b-a)}")

    # Two interleaved four-cycles with repeated vertices, plus isolated points.
    # Dense components trigger the automatic implicit route at the 2048-vertex probe.
    path = f"{TMP}/components.txt"
    with open(path, "w") as fh:
        for _ in range(96):
            for x, y in ((0, 0), (1, 0), (1, 1), (0, 1)):
                for offset in (0, 10):
                    print(x + offset, y, file=fh)
        for i in range(1280):
            print(20 + 2 * i, 0, file=fh)
    expected = Counter({(0, "0", "1"): 6, (0, "0", "inf"): 1282, (1, "1", "inf"): 2})
    for threads, prime in ((1, 2), (2, 3)):
        args = ["points", path, "--maxdim", "2", "--threshold", "1", "--threads", str(threads)]
        out, stats = run([PERSEUS] + args + ["--coeff", str(prime), "--print", "ripser", "--stats"])
        check(f"G disconnected Rips, {threads} workers F_{prime}", parse_ripser(out) == expected)
        check(f"G component split used with {threads} workers",
              f"2 nontrivial components, 1280 isolated vertices, {threads} workers" in stats, stats)
        check("G component field timing", "[solve " in stats and "field run " in stats, stats)

    # The barycentric subdivision of RP^2 is flag. Wedge a contractible clique at
    # one vertex of each copy to exercise the automatic split at this graph size.
    with open(os.path.join(HERE, "..", "examples", "data", "rp2_simplices.txt")) as fh:
        triangles = [line.split()[:-1] for line in fh if not line.startswith("#")]
    faces = sorted({frozenset(face) for t in triangles for k in (1, 2, 3)
                    for face in itertools.combinations(t, k)}, key=lambda f: (len(f), sorted(f)))
    path = f"{TMP}/component_rp2.flag"
    with open(path, "w") as fh:
        for v in range(2048):
            print(v, 0, file=fh)
        for offset in (0, 512):
            for i, face in enumerate(faces):
                for j, other in enumerate(faces[:i]):
                    if other < face:
                        weight = 0.5 if len(face) == 2 and face <= {"0", "1", "3"} else 1
                        print(offset + j, offset + i, weight, file=fh)
            for i, j in itertools.combinations([0] + list(range(len(faces), len(faces) + 260)), 2):
                print(offset + i, offset + j, 1, file=fh)
    args = ["flag", path, "--maxdim", "2", "--threshold", "1", "--print", "json"]
    reference = json.loads(run([PERSEUS] + args + ["--proved", "--engine", "implicit"])[0])
    for threads in (1, 2):
        out, stats = run([PERSEUS] + args + ["--threads", str(threads), "--stats"])
        got = json.loads(out)
        check(f"G proved RP^2 split used with {threads} workers",
              "2 nontrivial components," in stats and f"{threads} workers]" in stats, stats)
        check("G candidate rerun reuses each collapsed graph", stats.count("[edge collapse:") == 2, stats)
        check("G proved split matches unsplit certificate",
              all(got[k] == reference[k] for k in ("barcode", "bad_primes", "proof", "complete")))
        check("G proved split detects both copies of projective-plane torsion",
              got["proof"] and got["complete"] and got["candidates"] == [2]
              and [(b["p"], b["first_degree"], b["only_over_Fp"], b["only_over_Q"]) for b in got["bad_primes"]]
              == [(2, 1, [[], [[0.5, None]] * 2, [[1, None]] * 2], [[], [[0.5, 1]] * 2, []])])
    candidates, stats = run([PERSEUS] + args + ["--no-classify", "--stats"])
    cert = json.loads(candidates)
    check("G component candidate-only certificate", "nontrivial components," in stats and
          cert["proof"] and cert["candidates"] == [2] and cert["bad_primes_status"] == "candidates only"
          and not cert["complete"] and cert["bad_primes"] == [])


if __name__ == "__main__":
    which = sys.argv[1:] or ["A", "B", "C", "D", "E", "F", "G"]
    for s in which:
        print(f"=== section {s}", flush=True)
        if s.startswith("A:"): section_A(s[2:].split(","))
        else: globals()["section_" + s]()
    print(f"\n{results['pass']} passed, {results['fail']} failed")
    for f in failures[:40]: print("  ", f)
    sys.exit(1 if results["fail"] else 0)
