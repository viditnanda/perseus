#!/usr/bin/env python3
"""Cone-level engine vs the edge-collapse pipeline: GUDHI collapse_edges (Boissonnat-Pritam) then Ripser on
the collapsed flag filtration (sparse format). F_2 throughout; barcodes compared exactly."""
import numpy as np, gudhi, subprocess, time, sys, json, os
from collections import Counter
B = "/home/claude/bench/bin"; RIP = "/home/claude/ext/ripser/ripser"; P6 = "/home/claude/perseus6/perseus"
def parse(out):
    c = []; k = None
    for l in out.splitlines():
        if l.startswith("persistence intervals in dim"): k = int(l.split()[-1].rstrip(":")); continue
        l = l.strip()
        if k is not None and l.startswith("["):
            a, b = l[1:].rstrip(")").split(",")
            c.append((k, float(a), float("inf") if not b.strip() else float(b)))
    return sorted(c)
def same(x, y, tol=2e-5):
    # multiset match within a relative tolerance (Ripser prints 6 significant digits)
    import bisect
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
rows = []
for f, dim, thr in [("sphere", 2, None), ("o3_1024", 3, 1.8), ("o3_4096", 3, 1.4), ("random16", 7, None), ("fractal", 2, None), ("dragon", 1, None)]:
    raw = np.fromfile(f"{B}/{f}.bin", dtype=np.float32); n = int(round((1 + (1 + 8 * len(raw)) ** 0.5) / 2))
    D = np.zeros((n, n), dtype=np.float32); D[np.tril_indices(n, -1)] = raw; D = D + D.T
    if thr is None: thr = float(np.min(np.max(D, axis=1)))   # enclosing radius, as Ripser and Perseus use
    t0 = time.time()
    st = gudhi.RipsComplex(distance_matrix=D.astype(float), max_edge_length=thr).create_simplex_tree(max_dimension=1)
    t1 = time.time(); ne0 = st.num_simplices() - n
    for it in range(10):
        before = st.num_simplices(); st.collapse_edges(nb_iterations=1)
        if st.num_simplices() == before: break
    t2 = time.time(); ne1 = st.num_simplices() - n
    path = f"/tmp/coll_{f}.txt"
    with open(path, "w") as fh:
        for s, v in st.get_skeleton(1):
            if len(s) == 2: fh.write(f"{s[0]} {s[1]} {v}\n")
            else: fh.write(f"{s[0]} {s[0]} {v}\n")
    t3 = time.time()
    r = subprocess.run([RIP, "--format", "sparse", "--dim", str(dim), path], capture_output=True, text=True); t4 = time.time()
    t5 = time.time(); q = subprocess.run([P6, "bin32", f"{B}/{f}.bin", "--maxdim", str(dim), "--threshold", str(thr), "--coeff", "2"], capture_output=True, text=True); t6 = time.time()
    a, b = parse(r.stdout), parse(q.stdout)
    row = dict(case=f"{f} H<={dim}", edges=ne0, edges_after_collapse=ne1, build_s=round(t1 - t0, 2), collapse_s=round(t2 - t1, 2),
               ripser_on_collapsed_s=round(t4 - t3, 2), collapse_pipeline_total_s=round(t4 - t0 - (t3 - t2), 2), perseus6_s=round(t6 - t5, 2), agree=same(a, b))
    rows.append(row); print(json.dumps(row), flush=True)
json.dump(rows, open("/home/claude/perseus6/bench/vs_collapse.json", "w"), indent=1)
