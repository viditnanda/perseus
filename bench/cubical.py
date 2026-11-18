#!/usr/bin/env python3
"""Cubical (T-construction) race: Perseus 6 (discrete-gradient engine, and the classic implicit engine)
against CubicalRipser (tcripser.computePH) and GUDHI (CubicalPersistence), over F_2.

Every program runs in a fresh process started by a small wrapper (Linux keeps a process's peak RSS across
exec, so a process forked from a large parent would report the parent's memory); time is wall clock,
memory is the peak RSS of that process.
Perseus is timed from the command line (start, reading the .npy file, output suppressed); CubicalRipser
and GUDHI are timed around the in-process call, and their memory includes Python and the input array.
Barcodes are compared as multisets of positive-length intervals (Perseus vs CubicalRipser on every case,
Perseus engines against each other). The inputs are regenerated from a fixed seed.

usage: bench/cubical.py [out.json] [--quick] [--no-gudhi] [--cases a,b,...]
"""
import json, os, resource, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
P6 = os.environ.get("PERSEUS", os.path.join(HERE, "..", "perseus"))
DATA = os.environ.get("CUBDATA", "/tmp/cubical_bench")
os.makedirs(DATA, exist_ok=True)

def make(name):
    from scipy.ndimage import gaussian_filter
    kind, size = name.split("_"); n, d = map(int, size.split("x"))
    rng = np.random.default_rng(1000 * d + n + (0 if kind == "noise" else 7))
    A = rng.random((n,) * d)
    if kind == "smooth": A = gaussian_filter(A, {2: 4, 3: 2, 4: 1.5}[d])
    path = f"{DATA}/{name}.npy"
    if not os.path.exists(path): np.save(path, A)
    return path

CASES = ["smooth_1024x2", "noise_1024x2", "smooth_96x3", "noise_96x3", "smooth_24x4", "noise_24x4",
         "smooth_4096x2", "noise_4096x2", "smooth_192x3", "noise_192x3", "smooth_32x4", "noise_32x4", "smooth_48x4", "noise_48x4"]

WRAP = "import resource,subprocess,sys,time; t=time.time(); r=subprocess.run(sys.argv[1:], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True); print(time.time()-t, resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss, r.returncode)"

def perseus(path, engine, coeff="2"):
    out = subprocess.run([sys.executable, "-c", WRAP, P6, "npy", path, "--coeff", coeff, "--print", "none", "--engine", engine], capture_output=True, text=True).stdout.split()
    return float(out[0]), int(out[1]) / 1024, int(out[2])

def perseus_bars(path, engine):
    j = json.loads(subprocess.run([P6, "npy", path, "--coeff", "2", "--print", "json", "--engine", engine], capture_output=True, text=True).stdout)
    from collections import Counter
    return Counter((k, float(b), float("inf") if e is None else float(e)) for k, lst in enumerate(j["barcode"]) for b, e in lst if e is None or e > b)

INPROC = r'''
import sys, time, json, resource, numpy as np
from collections import Counter
which, path, res = sys.argv[1], sys.argv[2], sys.argv[3]
A = np.load(path); d = A.ndim
# imports before the clock starts (GUDHI's sklearn wrapper alone takes most of a second to import)
if which == "cripser": import tcripser
else: from gudhi.sklearn.cubical_persistence import CubicalPersistence
t = time.time()
if which == "cripser":
    r = tcripser.computePH(A, maxdim=d, location="none")
    dt = time.time() - t
    bars = Counter((int(x[0]), float(x[1]), float("inf") if x[2] >= 1.7e308 else float(x[2])) for x in r if x[2] > x[1])
else:
    dims = list(range(d))           # as in the README table: degrees 0..d-1 (H_d of a full grid is 0)
    out = CubicalPersistence(homology_dimensions=dims, homology_coeff_field=2, n_jobs=1).fit_transform([A])[0]
    dt = time.time() - t
    bars = Counter((k, float(b), float(e)) for k, arr in zip(dims, out) for b, e in arr if e > b)
if len(sys.argv) > 4: json.dump([[k, b, (None if e == float("inf") else e), c] for (k, b, e), c in bars.items()], open(sys.argv[4], "w"))
open(res, "w").write(str(dt))
'''

def other(which, path, dump=None):
    res = f"{DATA}/_timing.txt"
    if os.path.exists(res): os.remove(res)
    cmd = [sys.executable, "-c", WRAP, sys.executable, "-c", INPROC, which, path, res] + ([dump] if dump else [])
    w = subprocess.run(cmd, capture_output=True, text=True).stdout.split()
    if not os.path.exists(res) or int(w[2]) != 0: return None, None
    return float(open(res).read().split()[0]), int(w[1]) / 1024

def load_dump(fn):
    from collections import Counter
    return Counter({(k, b, float("inf") if e is None else e): c for k, b, e, c in json.load(open(fn))})

if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "--compare":
    path, dump = sys.argv[2], sys.argv[3]
    pm = perseus_bars(path, "morse"); pc = perseus_bars(path, "classic")
    out = {"morse_eq_classic": pm == pc, "intervals": sum(pm.values())}
    if os.path.exists(dump): out["morse_eq_cripser"] = pm == load_dump(dump)
    print(json.dumps(out)); sys.exit(0)

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    out = args[0] if args else os.path.join(HERE, "cubical.json")
    cases = CASES[:6] if "--quick" in sys.argv else CASES
    for a in sys.argv:
        if a.startswith("--cases="): cases = a.split("=", 1)[1].split(",")
    rows = []
    for name in cases:
        path = make(name)
        row = {"case": name}
        row["morse_s"], row["morse_MB"], _ = perseus(path, "morse")
        row["classic_s"], row["classic_MB"], _ = perseus(path, "classic")
        row["morse_Q_s"], row["morse_Q_MB"], _ = perseus(path, "morse", "Q")
        dump = f"{DATA}/{name}.cr.json"
        row["cripser_s"], row["cripser_MB"] = other("cripser", path, dump)
        if "--no-gudhi" not in sys.argv: row["gudhi_s"], row["gudhi_MB"] = other("gudhi", path)
        cmp = json.loads(subprocess.run([sys.executable, os.path.abspath(__file__), "--compare", path, dump], capture_output=True, text=True).stdout)
        row.update(cmp)
        rows.append(row); print(json.dumps(row), flush=True)
        json.dump(rows, open(out, "w"), indent=1)
