#!/usr/bin/env python3
"""Perseus 6 legacy formats vs the corresponding GUDHI constructions (and Perseus 5), over F_2.

For each case the same data are fed to
  * Perseus 6 CLI on the legacy text file (times include parsing),
  * Perseus 6 Python bindings on the in-memory data (same footing as GUDHI),
  * GUDHI (CubicalComplex / PeriodicCubicalComplex / SimplexTree / RipsComplex),
  * Perseus 5 on the legacy text file,
and the barcodes are compared as exact multisets (zero-length intervals dropped everywhere).
"""
import os, subprocess, sys, time, json, resource
from collections import Counter
from itertools import combinations
import numpy as np, gudhi, perseus

P6 = "/home/claude/perseus6/perseus"; P5 = "/home/claude/build/perseus"; TMP = "/tmp/vsg"; os.makedirs(TMP, exist_ok=True)
rng = np.random.default_rng(11)

def key(k, b, d): return (int(k), "%.6g" % b, "inf" if (d is None or np.isinf(d) or d == -1) else "%.6g" % d)
def from_gudhi(pers): return Counter(key(k, b, d) for k, (b, d) in pers if d > b)
def from_p6py(r): return Counter(key(k, b, d) for k, arr in enumerate(r["barcode"]) for b, d in arr)
def from_files(prefix, maxdim=8):
    c = Counter()
    for k in range(maxdim + 1):
        fn = f"{prefix}_{k}.txt"
        if os.path.exists(fn):
            for line in open(fn):
                t = line.split()
                if len(t) == 2: c[key(k, float(t[0]), float(t[1]))] += 1
    return c
def timed(f):
    t = time.time(); x = f(); return x, time.time() - t
def cli(prog, fmt, path, prefix):
    for fn in os.listdir(TMP):
        if fn.startswith(os.path.basename(prefix) + "_"): os.remove(os.path.join(TMP, fn))
    t = time.time(); subprocess.run([prog, fmt, path, prefix], stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=1800)
    return time.time() - t

rows = []
def report(name, size, bars, times):
    ref = bars["gudhi"]
    agree = {k: (v == ref) for k, v in bars.items() if k != "gudhi"}
    row = dict(case=name, size=size, intervals=sum(ref.values()), agree_with_gudhi=agree, seconds={k: round(v, 3) for k, v in times.items()})
    rows.append(row); print(json.dumps(row), flush=True)
    for k, v in bars.items():
        if k != "gudhi" and v != ref:
            print("   DIFF", k, "minus gudhi:", list((v - ref).items())[:5], " gudhi minus:", list((ref - v).items())[:5])

# ---------------------------------------------------------------- cubical (cubtop) and periodic cubical
def smooth_field(shape, nb=80):
    grids = np.meshgrid(*[np.arange(s) for s in shape], indexing="ij")
    f = np.zeros(shape)
    for _ in range(nb):
        c = [rng.random() * s for s in shape]; w = rng.random() * max(shape) / 12 + 2
        f += rng.random() * np.exp(-sum((g - ci) ** 2 for g, ci in zip(grids, c)) / (2 * w * w))
    f += 0.15 * rng.random(shape)
    return 1 + np.floor(254 * (f - f.min()) / (f.max() - f.min())).astype(int)

def cubical_case(name, shape, periodic):
    A = smooth_field(shape)
    # legacy cubtop lists top cubes with axis 0 fastest = Fortran order
    path = f"{TMP}/{name}.txt"
    ext = [(-s if p else s) for s, p in zip(shape, periodic)]
    with open(path, "w") as fh:
        fh.write("\n".join(map(str, [len(shape)] + ext)) + "\n"); np.savetxt(fh, A.flatten(order="F"), fmt="%d")
    bars, times = {}, {}
    def g():
        if any(periodic): cc = gudhi.PeriodicCubicalComplex(top_dimensional_cells=A.astype(float), periodic_dimensions=periodic)
        else:   # GUDHI's recommended entry point; in 2D it uses a dedicated image algorithm
            from gudhi.sklearn.cubical_persistence import CubicalPersistence
            dims = list(range(A.ndim))
            out = CubicalPersistence(homology_dimensions=dims, homology_coeff_field=2, n_jobs=1).fit_transform([A.astype(float)])[0]
            return [(k, (b, d)) for k, arr in zip(dims, out) for b, d in arr]
        return cc.persistence(homology_coeff_field=2)
    pg, times["gudhi"] = timed(g); bars["gudhi"] = from_gudhi(pg)
    r, times["perseus6_py"] = timed(lambda: perseus.cubical(A.astype(float), periodic=list(periodic), coeff=2)); bars["perseus6_py"] = from_p6py(r)
    times["perseus6_cli"] = cli(P6, "cubtop", path, f"{TMP}/{name}_p6"); bars["perseus6_cli"] = from_files(f"{TMP}/{name}_p6")
    if A.size <= 2_000_000:
        times["perseus5_cli"] = cli(P5, "cubtop", path, f"{TMP}/{name}_p5"); bars["perseus5_cli"] = from_files(f"{TMP}/{name}_p5")
    report(f"cubtop {name}" + (" periodic" if any(periodic) else ""), "x".join(map(str, shape)), bars, times)

# ---------------------------------------------------------------- simplicial (nmfsimtop): Linial-Meshulam process
def simtop_case(name, n, ntri):
    tris = list(combinations(range(n), 3)); idx = rng.permutation(len(tris))[:ntri]
    simp = [((a, b), 1) for a, b in combinations(range(n), 2)] + [(tris[i], 2 + t) for t, i in enumerate(idx)]
    path = f"{TMP}/{name}.txt"
    with open(path, "w") as fh:
        fh.write("1\n")
        for s, b in simp: fh.write(f"{len(s) - 1} " + " ".join(str(v + 1) for v in s) + f" {b}\n")
    bars, times = {}, {}
    def g():
        st = gudhi.SimplexTree()
        for s, b in simp: st.insert(list(s), filtration=b)
        return st.persistence(homology_coeff_field=2, persistence_dim_max=True)
    pg, times["gudhi"] = timed(g); bars["gudhi"] = from_gudhi(pg)
    r, times["perseus6_py"] = timed(lambda: perseus.simplices([(list(s), float(b)) for s, b in simp], coeff=2)); bars["perseus6_py"] = from_p6py(r)
    times["perseus6_cli"] = cli(P6, "nmfsimtop", path, f"{TMP}/{name}_p6"); bars["perseus6_cli"] = from_files(f"{TMP}/{name}_p6")
    times["perseus5_cli"] = cli(P5, "nmfsimtop", path, f"{TMP}/{name}_p5"); bars["perseus5_cli"] = from_files(f"{TMP}/{name}_p5")
    report(f"nmfsimtop {name}", f"n={n}, {len(simp)} listed simplices", bars, times)

# ---------------------------------------------------------------- distmat: Rips on step births
def distmat_case(name, n, g0, s, N, d):
    M = rng.random((n, n)) * 10; M = np.round((M + M.T) / 2, 3); np.fill_diagonal(M, 0)
    path = f"{TMP}/{name}.txt"
    with open(path, "w") as fh:
        fh.write(f"{n}\n{g0} {s} {N} {d}\n"); np.savetxt(fh, M, fmt="%.3f")
    steps = np.where(M < g0, 0, np.ceil((M - g0) / s)); steps[steps > N] = np.inf; np.fill_diagonal(steps, 0)
    bars, times = {}, {}
    def g():
        D = np.where(np.isinf(steps), 1e300, steps)
        st = gudhi.RipsComplex(distance_matrix=D, max_edge_length=float(N)).create_simplex_tree(max_dimension=d)
        return st.persistence(homology_coeff_field=2, persistence_dim_max=True)
    pg, times["gudhi"] = timed(g); bars["gudhi"] = from_gudhi(pg)
    times["perseus6_cli"] = cli(P6, "distmat", path, f"{TMP}/{name}_p6"); bars["perseus6_cli"] = from_files(f"{TMP}/{name}_p6")
    times["perseus5_cli"] = cli(P5, "distmat", path, f"{TMP}/{name}_p5"); bars["perseus5_cli"] = from_files(f"{TMP}/{name}_p5")
    report(f"distmat {name}", f"n={n}, simplices to dim {d}", bars, times)

# ---------------------------------------------------------------- brips: growing balls, vertices from step 1
def brips_case(name, n, g0, s, N):
    X = np.round(rng.random((n, 2)) * 4, 3); R = np.round(rng.random(n) * 0.3, 2)
    path = f"{TMP}/{name}.txt"
    with open(path, "w") as fh:
        fh.write(f"2\n{g0} {s} {N}\n")
        for x, r in zip(X, R): fh.write(f"{x[0]} {x[1]} {r}\n")
    D = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))
    steps = np.maximum(1, np.maximum(0, np.ceil(0.5 / s * (D / g0 - (R[:, None] + R[None]))))); np.fill_diagonal(steps, 0)
    bars, times = {}, {}
    def g():
        st = gudhi.RipsComplex(distance_matrix=steps, max_edge_length=float(N)).create_simplex_tree(max_dimension=2)
        for v in range(n): st.assign_filtration([v], 1.0)
        st.make_filtration_non_decreasing()
        return st.persistence(homology_coeff_field=2, persistence_dim_max=True)
    pg, times["gudhi"] = timed(g); bars["gudhi"] = from_gudhi(pg)
    times["perseus6_cli"] = cli(P6, "brips", path, f"{TMP}/{name}_p6"); bars["perseus6_cli"] = from_files(f"{TMP}/{name}_p6")
    times["perseus5_cli"] = cli(P5, "brips", path, f"{TMP}/{name}_p5"); bars["perseus5_cli"] = from_files(f"{TMP}/{name}_p5")
    report(f"brips {name}", f"n={n} points in R^2", bars, times)

if __name__ == "__main__":
    cubical_case("img256", (256, 256), (False, False))
    cubical_case("img1000", (1000, 1000), (False, False))
    cubical_case("img2000", (2000, 2000), (False, False))
    cubical_case("vol64", (64, 64, 64), (False, False, False))
    cubical_case("vol128", (128, 128, 128), (False, False, False))
    cubical_case("torus512", (512, 512), (True, True))
    cubical_case("cyl64", (64, 64, 64), (True, False, True))
    simtop_case("lm40", 40, 9880)
    simtop_case("lm70", 70, 30000)
    distmat_case("dm60", 60, 0.5, 0.25, 30, 3)
    distmat_case("dm120", 120, 0.5, 0.25, 12, 3)
    brips_case("br150", 150, 1.0, 0.1, 15)
    json.dump(rows, open(sys.argv[1] if len(sys.argv) > 1 else "vs_gudhi.json", "w"), indent=1)
