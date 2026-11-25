#!/usr/bin/env python3
"""Vietoris-Rips race, 25 Nov 2026: Perseus (Q certified and F_2, command line), Ripser 1.2 (command line),
RedZeD (redzed-tda 0.1.1, in-process) and GUDHI 3.13 edge collapse followed by Ripser, on Ripser's benchmark
inputs as float32 lower triangles in $BENCH (sphere_3_192, o3_1024, o3_4096, random_point_cloud_50_16 as
random16, the 512-point fractal with random edge weights from the PH roadmap as fractal, dragon_vrip 2000
points, and 180 points with uniformly random distances, seed 180). Barcodes are compared with Perseus's F_2
barcode to 2e-5 relative (Ripser prints 6 significant digits). Times: best of three below 10 s.

usage: PERSEUS=... RIPSER=... BENCH=... python3 bench/rips_race.py [out.json]
"""
import json, os, subprocess, sys, time, bisect
import numpy as np, gudhi
P6 = os.environ.get("PERSEUS", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "perseus"))
RIP = os.environ.get("RIPSER", "ripser"); B = os.environ.get("BENCH", "/tmp/rips_bench")
cases = [("sphere", 2, None), ("o3_1024", 3, 1.8), ("o3_4096", 1, 1.4), ("o3_4096", 3, 1.4), ("random16", 1, None),
         ("random16", 7, None), ("fractal", 1, None), ("fractal", 2, None), ("dragon", 1, None), ("randdist180", 1, None)]
def load(f):
    raw = np.fromfile(f"{B}/{f}.bin", dtype=np.float32); n = int(round((1 + (1 + 8 * len(raw)) ** 0.5) / 2))
    D = np.zeros((n, n), dtype=np.float64); D[np.tril_indices(n, -1)] = raw; return D + D.T
def timed(cmd, reps):
    best = None
    for _ in range(reps):
        p = subprocess.run(["/usr/bin/time", "-f", "%e %M"] + cmd, capture_output=True, text=True, timeout=1800)
        t, m = p.stderr.strip().splitlines()[-1].split(); t = float(t); m = int(m) / 1024
        if best is None or t < best[0]: best = (t, m, p.stdout)
        if t > 10: break
    return best
def parse_ripser(out):
    c = []; k = None
    for l in out.splitlines():
        if l.startswith("persistence intervals in dim"): k = int(l.split()[-1].rstrip(":")); continue
        l = l.strip()
        if k is not None and l.startswith("["):
            a, b = l[1:].rstrip(")").split(","); c.append((k, float(a), float("inf") if not b.strip() else float(b)))
    return [x for x in c if x[1] != x[2]]
def same(x, y, tol=2e-5):
    if len(x) != len(y): return False
    close = lambda p, q: p == q or (p != float("inf") and q != float("inf") and abs(p - q) <= tol * max(1.0, abs(p)))
    for k in set(t[0] for t in x) | set(t[0] for t in y):
        X = sorted((b, d) for kk, b, d in x if kk == k); Y = sorted((b, d) for kk, b, d in y if kk == k)
        if len(X) != len(Y): return False
        births = [b for b, _ in Y]; used = [False] * len(Y)
        for b, d in X:
            i = bisect.bisect_left(births, b - tol * max(1.0, abs(b)) - 1e-12); found = False
            while i < len(Y) and births[i] <= b + tol * max(1.0, abs(b)) + 1e-12:
                if not used[i] and close(Y[i][1], d): used[i] = True; found = True; break
                i += 1
            if not found: return False
    return True
RZ = r'''
import sys, time, resource, json, numpy as np, redzed_tda as rz
raw = np.fromfile(sys.argv[1], dtype=np.float32).astype(np.float64)
n = int(round((1 + (1 + 8 * len(raw)) ** 0.5) / 2)); D = np.zeros((n, n)); D[np.tril_indices(n, -1)] = raw; D = D + D.T
dim = int(sys.argv[2]); thr = float(sys.argv[3])
t = time.perf_counter(); res = rz.rips(D, maxdim=dim, threshold=thr); t = time.perf_counter() - t
m = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
bars = [[k, float(b), float(d)] for k, A in enumerate(res) for b, d in np.asarray(A).reshape(-1, 2) if b != d]
print(json.dumps({"t": t, "m": m, "bars": bars}))
'''
open("/tmp/rz_run.py", "w").write(RZ)
rows = []
for f, dim, thr in cases:
    D = load(f); n = len(D)
    if thr is None: thr = float(np.min(np.max(D, axis=1)))
    T = ["--threshold", repr(thr)]
    path = f"{B}/{f}.bin"
    q = timed([P6, "bin32", path, "--maxdim", str(dim), "--print", "none"] + T, 3)
    f2 = timed([P6, "bin32", path, "--maxdim", str(dim), "--coeff", "2"] + T, 3)
    rp = timed([RIP, "--format", "binary", "--dim", str(dim), "--threshold", repr(thr), path], 3)
    rzb = None
    for _ in range(3):
        p = subprocess.run([sys.executable, "/tmp/rz_run.py", path, str(dim), repr(thr)], capture_output=True, text=True, timeout=1800)
        if p.returncode: r = {"t": float("nan"), "m": float("nan"), "bars": [], "error": p.stderr[-300:]}; rzb = r; break
        r = json.loads(p.stdout)
        if rzb is None or r["t"] < rzb["t"]: rzb = r
        if r["t"] > 10: break
    a = parse_ripser(f2[2]); b = parse_ripser(rp[2]); c = [tuple(x) for x in rzb["bars"]]
    row = dict(case=f"{f} H<={dim}", n=n, threshold=thr, perseus_Q_s=round(q[0], 2), perseus_Q_MB=round(q[1]),
               perseus_F2_s=round(f2[0], 2), perseus_F2_MB=round(f2[1]), ripser_s=round(rp[0], 2), ripser_MB=round(rp[1]),
               redzed_s=round(rzb["t"], 2), redzed_MB=rzb["m"] if rzb["m"] != rzb["m"] else round(rzb["m"]), agree_ripser=same(a, b), agree_redzed=same(a, c), redzed_error=rzb.get("error"))
    if dim >= 2 or f == "dragon":   # edge collapse (GUDHI) + Ripser on the collapsed graph
        t0 = time.time()
        st = gudhi.RipsComplex(distance_matrix=D, max_edge_length=thr).create_simplex_tree(max_dimension=1)
        e0 = st.num_simplices() - n
        for it in range(10):
            before = st.num_simplices(); st.collapse_edges(nb_iterations=1)
            if st.num_simplices() == before: break
        e1 = st.num_simplices() - n
        path2 = f"/tmp/coll_{f}.txt"
        with open(path2, "w") as fh:
            for s, v in st.get_skeleton(1):
                fh.write(f"{s[0]} {s[1]} {v}\n" if len(s) == 2 else f"{s[0]} {s[0]} {v}\n")
        tw = time.time()
        r2 = subprocess.run([RIP, "--format", "sparse", "--dim", str(dim), path2], capture_output=True, text=True); t1 = time.time()
        row.update(edges=e0, edges_after_collapse=e1, collapse_ripser_s=round((tw - t0) + (t1 - tw) - 0, 2), agree_collapse=same(a, parse_ripser(r2.stdout)))
    rows.append(row); print(json.dumps(row), flush=True)
json.dump(rows, open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "rips_race.json"), "w"), indent=1)
