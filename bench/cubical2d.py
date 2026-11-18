#!/usr/bin/env python3
"""2D race: Perseus (planar engine, and the gradient engine for reference) against GUDHI and
CubicalRipser on top-cell (T-construction) filtrations of 2D images, over F_2.

Like with like: every tool is timed around its in-process Python call, after its module is imported (input already in memory as a
NumPy array; output returned as arrays), each in a fresh process started by a small wrapper. Memory is the
process's peak RSS read right after the call (so it includes the interpreter, the modules and the input
array, but not the comparison of barcodes done afterwards); Linux keeps peak RSS across exec, so the process
must not be forked from a large parent.
Perseus's command line (process start, reading the .npy file, output suppressed) is timed separately.
GUDHI's CubicalPersistence is asked for degrees 0 and 1 (its 2D image algorithm); CubicalRipser
(tcripser.computePH) and Perseus compute all degrees. Barcodes (positive-length intervals) are compared.

usage: bench/cubical2d.py [out.json] [--cases=a,b,...] [--repeat=n] [--tools=planar,gudhi,...]
  (at 8192^2 the gradient engine and CubicalRipser are skipped: they need more memory than the 7 GB test machine)
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
import sys, time, json, resource, numpy as np
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
elif which == "planar_Q":
    r = perseus.cubical(A)                                   # Q, certified
elif which == "cripser":
    out = tcripser.computePH(A, maxdim=d, location="none")
else:
    out = CubicalPersistence(homology_dimensions=[0, 1], homology_coeff_field=2, n_jobs=1).fit_transform([A])[0]
dt = time.perf_counter() - t; peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
open(res, "w").write(f"{dt} {peak}")
if dump != "-":       # positive-length intervals as a sorted (degree, birth, death) array
    if which.startswith("planar") or which == "morse": parts = [(k, np.asarray(a, dtype=float).reshape(-1, 2)) for k, a in enumerate(r["barcode"])]
    elif which == "cripser":
        out = np.asarray(out); dth = np.where(out[:, 2] >= 1.7e308, np.inf, out[:, 2])
        parts = [(k, np.stack([out[out[:, 0] == k, 1], dth[out[:, 0] == k]], axis=1)) for k in range(d + 1)]
    else: parts = [(k, np.asarray(a, dtype=float).reshape(-1, 2)) for k, a in zip([0, 1], out)]
    B = np.concatenate([np.column_stack([np.full(len(a), k, dtype=float), a]) for k, a in parts if len(a)] or [np.zeros((0, 3))])
    B = B[B[:, 2] > B[:, 1]]
    B = B[np.lexsort((B[:, 2], B[:, 1], B[:, 0]))]
    np.save(dump, B)
if which == "planar_Q": open(res + ".cert", "w").write(json.dumps(r.get("certificate", {}), default=str))
'''

def inproc(which, path, dump="-"):
    res = f"{DATA}/_t2d.txt"
    if os.path.exists(res): os.remove(res)
    w = subprocess.run([sys.executable, "-c", WRAP, sys.executable, "-c", INPROC, which, path, res, dump, PYPATH], capture_output=True, text=True)
    parts = w.stdout.split()
    if not os.path.exists(res) or len(parts) < 2 or int(parts[1]) != 0: return None, None, w.stderr[-300:]
    dt, peak = open(res).read().split()
    return float(dt), int(peak) / 1024, None

def cli(path):
    w = subprocess.run([sys.executable, "-c", WRAP, "/usr/bin/env", "python3", "-c",
                        f"import subprocess,time,sys; t=time.perf_counter(); subprocess.run(['{P6}','npy','{path}','--coeff','2','--print','none'],check=True); open('{DATA}/_cli.txt','w').write(str(time.perf_counter()-t))"],
                       capture_output=True, text=True)
    return float(open(f"{DATA}/_cli.txt").read()), int(w.stdout.split()[0]) / 1024

def load(fn):
    return np.load(fn)

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    out = args[0] if args else os.path.join(HERE, "cubical2d.json")
    cases = CASES; rep = 1; tools = ["planar", "planar_Q", "morse", "gudhi", "cripser"]
    for a in sys.argv:
        if a.startswith("--cases="): cases = a.split("=", 1)[1].split(",")
        if a.startswith("--repeat="): rep = int(a.split("=", 1)[1])
        if a.startswith("--tools="): tools = a.split("=", 1)[1].split(",")
    rows = []
    for name in cases:
        path = make(name); row = {"case": name}
        dumps = {}
        if os.path.exists(f"{DATA}/_t2d.txt.cert"): os.remove(f"{DATA}/_t2d.txt.cert")
        for tool in tools:
            if name.split("_")[1].startswith("8192") and tool in ("morse", "cripser"):
                row[tool + "_skipped"] = "memory: beyond this machine's 7 GB at 8192^2"; continue
            best = None
            for r in range(rep):
                dump = f"{DATA}/{name}.{tool}.bars.npy" if r == 0 and tool in ("planar", "gudhi", "cripser", "morse") else "-"
                t, mb, err = inproc(tool, path, dump)
                if t is None: row[tool + "_error"] = err; break
                if best is None or t < best[0]: best = (t, mb)
                if dump != "-": dumps[tool] = dump
            if best: row[tool + "_s"], row[tool + "_MB"] = best
            print(name, tool, best, flush=True)
        try: row["cli_s"], row["cli_MB"] = cli(path)
        except Exception as e: row["cli_error"] = str(e)
        if "planar" in dumps:
            P = load(dumps["planar"]); row["intervals"] = len(P)
            for tool in ("morse", "gudhi", "cripser"):
                if tool in dumps: row["planar_eq_" + tool] = bool(np.array_equal(P, load(dumps[tool])))   # a 2D grid has no H2 bars
        if os.path.exists(f"{DATA}/_t2d.txt.cert"): row["certificate"] = json.loads(open(f"{DATA}/_t2d.txt.cert").read())
        rows.append(row); print(json.dumps(row), flush=True)
        json.dump(rows, open(out, "w"), indent=1)
