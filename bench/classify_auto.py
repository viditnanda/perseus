#!/usr/bin/env python3
"""How the small candidate primes are classified on explicit complexes: --classify-by rerun (one rerun per
prime), sweep (batched Howell sweeps, torsion.hpp) and auto (a timed sample of reruns sets a budget; the
sweeps go ahead only if their reference sweeps, run under that budget, predict they cost less). auto must be
within noise of the better of the two on every instance, and all three must find identical bad primes.

usage: python3 bench/classify_auto.py [out.json]     (writes inputs to $PERSEUS_BENCH_TMP or /tmp/perseus_classify)
"""
import json, os, re, subprocess, sys, time, itertools
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PERSEUS = os.environ.get("PERSEUS", os.path.join(ROOT, "perseus"))
T = os.environ.get("PERSEUS_BENCH_TMP", "/tmp/perseus_classify"); os.makedirs(T, exist_ok=True)
THREADS = os.environ.get("THREADS", "2")


def lm(n, m, seed, d=2, distinct=False):
    path = f"{T}/lm{'d' if distinct else ''}{d}_{n}_{m}_{seed}.txt"
    if os.path.exists(path): return path
    rng = np.random.default_rng(seed)
    with open(path, "w") as f:
        if not distinct:
            for s in itertools.combinations(range(n), d): f.write(" ".join(map(str, s)) + " 1\n")
            tops = list(itertools.combinations(range(n), d + 1)); idx = rng.permutation(len(tops))[:m]
            for i, j in enumerate(idx): f.write(" ".join(map(str, tops[j])) + f" {i + 2}\n")
        else:
            E = {e: rng.random() for e in itertools.combinations(range(n), 2)}
            for (a, b), w in E.items(): f.write(f"{a} {b} {w:.12f}\n")
            tris = list(itertools.combinations(range(n), 3)); idx = rng.permutation(len(tris))[:m]
            for j in idx:
                a, b, c = tris[j]; me = max(E[(a, b)], E[(a, c)], E[(b, c)])
                f.write(f"{a} {b} {c} {max(rng.random() * 1.3, me + 1e-7 * (1 + rng.random())):.12f}\n")
    return path


def rp2_points(n, seed):
    path = f"{T}/rp2_{n}_{seed}.txt"
    if not os.path.exists(path):
        rng = np.random.default_rng(seed)
        v = rng.normal(size=(n, 3)); v /= np.linalg.norm(v, axis=1)[:, None]; x, y, z = v.T
        np.savetxt(path, np.c_[x * x - y * y, 2 * x * y, 2 * x * z, 2 * y * z], fmt="%.10f")
    return path


def mixed(npts, r, seed, ln=40, lmm=700):
    """the Rips complex of npts points on RP^2 up to dimension 3 at scale r, plus a disjoint Linial-Meshulam
    component with distinct interleaved values (its candidates)"""
    path = f"{T}/mixed_{npts}_{r}_{seed}" + ("" if ln == 40 else f"_lm{ln}") + ".txt"
    if os.path.exists(path): return path
    X = np.loadtxt(rp2_points(npts, seed))
    D = np.linalg.norm(X[:, None] - X[None], axis=2)
    adj = [set(np.nonzero((D[i] <= r) & (np.arange(npts) > i))[0]) for i in range(npts)]
    with open(path, "w") as f:
        def emit(s):
            f.write(" ".join(map(str, s)) + f" {max(D[a, b] for a, b in itertools.combinations(s, 2)):.10f}\n")
        for i in range(npts):
            for j in adj[i]:
                emit((i, j))
                for k in adj[i] & adj[j]:
                    emit((i, j, k))
                    for l in adj[i] & adj[j] & adj[k]: emit((i, j, k, l))
        off = 10 * npts
        for line in open(lm(ln, lmm, ln, distinct=True)):
            t = line.split(); f.write(" ".join(str(int(x) + off) for x in t[:-1]) + f" {t[-1]}\n")
    return path


def chain(n1, n2, seed):
    """a filtered chain complex over Z, 2-cells with random integer boundaries: many candidates"""
    path = f"{T}/chain_{n1}_{n2}_{seed}.cells"
    if os.path.exists(path): return path
    rng = np.random.default_rng(seed)
    coefs = [1, -1, 2, -2, 3, 5, 7, 11, 13, 6, 10, 9]
    with open(path, "w") as f:
        b1 = rng.random(n1)
        for b in b1: f.write(f"1 {b:.10f}\n")
        for _ in range(n2):
            fs = rng.choice(n1, size=rng.integers(2, 5), replace=False)
            v = max(b1[fs]) + rng.random()
            f.write(f"2 {v:.10f} " + " ".join(f"{a} {coefs[rng.integers(len(coefs))]}" for a in fs) + "\n")
    return path


def cubical(n, seed):
    path = f"{T}/noise_{n}_{seed}.npy"
    if not os.path.exists(path): np.save(path, np.random.default_rng(seed).random((n, n, n)))
    return path


CASES = [
    ("LM n=70, tied", ["simplices", lambda: lm(70, 2400, 70), "--maxdim", "2"]),
    ("LM n=90, tied", ["simplices", lambda: lm(90, 4000, 90), "--maxdim", "2"]),
    ("LM n=40, distinct values", ["simplices", lambda: lm(40, 700, 40, distinct=True), "--maxdim", "2"]),
    ("LM n=70, distinct values", ["simplices", lambda: lm(70, 2400, 70, distinct=True), "--maxdim", "2"]),
    ("LM d=3 n=25, tied", ["simplices", lambda: lm(25, 2700, 25, d=3), "--maxdim", "3"]),
    ("random chain complex 400/900", ["cells", lambda: chain(400, 900, 1), "--maxdim", "2"]),
    ("RP^2 Rips 300 points, sparse route, H<=2", ["points", lambda: rp2_points(300, 4), "--maxdim", "2", "--engine", "sparse"]),
    ("RP^2 Rips 300 points (dim <= 3, 0.9M tetrahedra) + LM n=40 component", ["simplices", lambda: mixed(300, 1.0, 5), "--maxdim", "2"]),
    ("RP^2 Rips 300 points (dim <= 3, 0.9M tetrahedra) + LM n=70 component", ["simplices", lambda: mixed(300, 1.0, 5, 70, 2400), "--maxdim", "2"]),
    ("cubical noise 40^3, explicit engine", ["npy", lambda: cubical(40, 1), "--maxdim", "2", "--engine", "explicit"]),
    ("fixture crt_alias_1", ["cells", lambda: os.path.join(ROOT, "tests/fixtures/crt_alias_1.cells")]),
    ("fixture crt_group", ["cells", lambda: os.path.join(ROOT, "tests/fixtures/crt_group.cells")]),
]


def run(args):
    env = dict(os.environ, PERSEUS_DEBUG="1")
    t0 = time.perf_counter()
    p = subprocess.run([PERSEUS] + args + ["--print", "json", "--threads", THREADS], capture_output=True, text=True, env=env)
    wall = time.perf_counter() - t0
    if p.returncode: raise RuntimeError(p.stderr[-2000:])
    j = json.loads(p.stdout)
    m = re.search(r"small candidates ([0-9.]+) s wall", p.stderr)
    sw = re.search(r"sweep classification (used|not used): (.*)\]", p.stderr)
    return wall, j, (float(m[1]) if m else 0.0), (sw[1] + ": " + sw[2] if sw else "no sweep attempted")


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else None
    rows = []
    for name, spec in CASES:
        args = [spec[0], spec[1]()] + spec[2:]
        row = {"case": name}
        res = {}
        for mode in ("rerun", "sweep", "auto"):
            wall, j, small, note = run(args + ["--classify-by", mode])
            res[mode] = j
            row[mode] = {"total_s": round(wall, 3), "small_s": small, "note": note}
        row["candidates"] = len(res["auto"]["candidates"])
        row["bad"] = len(res["auto"]["bad_primes"])
        row["identical"] = res["rerun"]["bad_primes"] == res["sweep"]["bad_primes"] == res["auto"]["bad_primes"]
        best = min(row["rerun"]["small_s"], row["sweep"]["small_s"])
        row["auto_vs_best_small"] = round(row["auto"]["small_s"] - best, 3)
        rows.append(row)
        print(f"{name}: cands {row['candidates']}, bad {row['bad']}, identical {row['identical']} | "
              f"small: rerun {row['rerun']['small_s']:.2f} s, sweep {row['sweep']['small_s']:.2f} s, auto {row['auto']['small_s']:.2f} s "
              f"| total: {row['rerun']['total_s']:.2f} / {row['sweep']['total_s']:.2f} / {row['auto']['total_s']:.2f} s | auto: {row['auto']['note']}",
              flush=True)
    if out:
        with open(out, "w") as f: json.dump(rows, f, indent=1)


if __name__ == "__main__":
    main()
