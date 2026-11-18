#!/usr/bin/env python3
"""2D race: Perseus (planar engine, and the gradient engine for reference) against GUDHI and
CubicalRipser on top-cell (T-construction) filtrations of 2D images, over F_2.

Like with like: every tool is timed around its in-process Python call, after its module is imported (input already in memory as a
NumPy array; output returned as arrays), each in a fresh process started by a small wrapper, whose peak
RSS is reported (Linux keeps peak RSS across exec, so the process must not be forked from a large parent).
Perseus's command line (process start, reading the .npy file, output suppressed) is timed separately.
GUDHI's CubicalPersistence is asked for degrees 0 and 1 (its 2D image algorithm); CubicalRipser
(tcripser.computePH) and Perseus compute all degrees. Barcodes (positive-length intervals) are compared.

usage: bench/cubical2d.py [out.json] [--cases a,b,...] [--repeat n]
  PERSEUS_PY: directory to put first on sys.path for the perseus module (default: the installed one)
"""
import json, os, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
P6 = os.environ.get("PERSEUS", os.path.join(HERE, "..", "perseus"))
DATA = os.environ.get("CUBDATA", "/tmp/cubical_bench")
PYPATH = os.environ.get("PERSEUS_PY", "")
os.makedirs(DATA, exist_ok=True)

def make(name):
    path = f"{DATA}/{name}.npy"
    if os.path.exists(path): return path
    if name == "photo_512x2":
        import skimage.data
        A = skimage.data.camera().astype(np.float64)          # 512 x 512, 8-bit grey levels: heavy ties
    else:
        from scipy.ndimage import gaussian_filter
        kind, size = name.split("_"); n = int(size.split("x")[0])
        rng = np.random.default_rng(2000 + n + (0 if kind == "noise" else 7))
        A = rng.random((n, n))
        if kind == "smooth": A = gaussian_filter(A, 4)
    np.save(path, A)
    return path

CASES = ["photo_512x2", "smooth_1024x2", "noise_1024x2", "smooth_4096x2", "noise_4096x2", "smooth_8192x2", "noise_8192x2"]

WRAP = "import resource,subprocess,sys; r=subprocess.run(sys.argv[1:], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True); print(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss, r.returncode); sys.stderr.write(r.stderr[-2000:])"

INPROC = r'''
import sys, time, json, numpy as np
from collections import Counter
which, path, res, dump, pypath = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
if pypath: sys.path.insert(0, pypath)
A = np.load(path); d = A.ndim
# imports before the clock starts (GUDHI's sklearn wrapper alone takes most of a second to import)
if which.startswith("planar") or which == "morse": import perseus
elif which == "cripser": import tcripser
else: from gudhi.sklearn.cubical_persistence import CubicalPersistence
t = time.perf_counter()
if which in ("planar", "morse"):
    r = perseus.cubical(A, coeff=2, engine=which)
    dt = time.perf_counter() - t
    bars = Counter((k, float(b), float(e)) for k, arr in enumerate(r["barcode"]) for b, e in np.asarray(arr).reshape(-1, 2) if e > b)
elif which == "planar_Q":
    r = perseus.cubical(A)                                   # Q, certified
    dt = time.perf_counter() - t
    bars = Counter((k, float(b), float(e)) for k, arr in enumerate(r["barcode"]) for b, e in np.asarray(arr).reshape(-1, 2) if e > b)
    open(res + ".cert", "w").write(json.dumps(r.get("certificate", {}), default=str))
elif which == "cripser":
    out = tcripser.computePH(A, maxdim=d, location="none")
    dt = time.perf_counter() - t
    bars = Counter((int(x[0]), float(x[1]), float("inf") if x[2] >= 1.7e308 else float(x[2])) for x in out if x[2] > x[1])
else:
    out = CubicalPersistence(homology_dimensions=[0, 1], homology_coeff_field=2, n_jobs=1).fit_transform([A])[0]
    dt = time.perf_counter() - t
    bars = Counter((k, float(b), float(e)) for k, arr in zip([0, 1], out) for b, e in arr if e > b)
if dump != "-": json.dump([[k, b, (None if e == float("inf") else e), c] for (k, b, e), c in bars.items()], open(dump, "w"))
open(res, "w").write(str(dt))
'''

def inproc(which, path, dump="-"):
    res = f"{DATA}/_t2d.txt"
    if os.path.exists(res): os.remove(res)
    w = subprocess.run([sys.executable, "-c", WRAP, sys.executable, "-c", INPROC, which, path, res, dump, PYPATH], capture_output=True, text=True)
    parts = w.stdout.split()
    if not os.path.exists(res) or len(parts) < 2 or int(parts[1]) != 0: return None, None, w.stderr[-300:]
    return float(open(res).read()), int(parts[0]) / 1024, None

def cli(path):
    w = subprocess.run([sys.executable, "-c", WRAP, "/usr/bin/env", "python3", "-c",
                        f"import subprocess,time,sys; t=time.perf_counter(); subprocess.run(['{P6}','npy','{path}','--coeff','2','--print','none'],check=True); open('{DATA}/_cli.txt','w').write(str(time.perf_counter()-t))"],
                       capture_output=True, text=True)
    return float(open(f"{DATA}/_cli.txt").read()), int(w.stdout.split()[0]) / 1024

def load(fn):
    from collections import Counter
    return Counter({(k, b, float("inf") if e is None else e): c for k, b, e, c in json.load(open(fn))})

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    out = args[0] if args else os.path.join(HERE, "cubical2d.json")
    cases = CASES; rep = 1
    for a in sys.argv:
        if a.startswith("--cases="): cases = a.split("=", 1)[1].split(",")
        if a.startswith("--repeat="): rep = int(a.split("=", 1)[1])
    rows = []
    for name in cases:
        path = make(name); row = {"case": name}
        dumps = {}
        if os.path.exists(f"{DATA}/_t2d.txt.cert"): os.remove(f"{DATA}/_t2d.txt.cert")
        for tool in ["planar", "planar_Q", "morse", "gudhi", "cripser"]:
            best = None
            for r in range(rep):
                dump = f"{DATA}/{name}.{tool}.json" if r == 0 and tool in ("planar", "gudhi", "cripser", "morse") else "-"
                t, mb, err = inproc(tool, path, dump)
                if t is None: row[tool + "_error"] = err; break
                if best is None or t < best[0]: best = (t, mb)
                if dump != "-": dumps[tool] = dump
            if best: row[tool + "_s"], row[tool + "_MB"] = best
            print(name, tool, best, flush=True)
        try: row["cli_s"], row["cli_MB"] = cli(path)
        except Exception as e: row["cli_error"] = str(e)
        if "planar" in dumps:
            P = load(dumps["planar"]); row["intervals"] = sum(P.values())
            for tool in ("morse", "gudhi", "cripser"):
                if tool in dumps: row["planar_eq_" + tool] = P == load(dumps[tool])   # a 2D grid has no H2 bars
        if os.path.exists(f"{DATA}/_t2d.txt.cert"): row["certificate"] = json.loads(open(f"{DATA}/_t2d.txt.cert").read())
        rows.append(row); print(json.dumps(row), flush=True)
        json.dump(rows, open(out, "w"), indent=1)
