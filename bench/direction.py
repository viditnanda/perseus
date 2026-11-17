#!/usr/bin/env python3
"""Which reduction direction suits which explicit complex? Times (F_2) of cohomology / homology / forward."""
import subprocess, re, sys, os, random, json
from itertools import combinations
import numpy as np
P6 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "perseus"); T = "/tmp/dir"
rng = random.Random(4)

def lm(name, n, d, frac):
    faces = list(combinations(range(1, n + 1), d + 1)); rng.shuffle(faces); faces = faces[:int(frac * len(faces))]
    with open(f"{T}/{name}.txt", "w") as f:
        f.write("1\n")
        for s in combinations(range(1, n + 1), d): f.write(f"{d - 1} " + " ".join(map(str, s)) + " 1\n")
        for i, s in enumerate(faces): f.write(f"{d} " + " ".join(map(str, s)) + f" {i + 2}\n")
    return ("nmfsimtop", f"{T}/{name}.txt", [])
def randsimp(name, n, nf, dmax):
    with open(f"{T}/{name}.txt", "w") as f:
        f.write("1\n")
        for _ in range(nf):
            d = rng.randint(1, dmax); s = rng.sample(range(1, n + 1), d + 1)
            f.write(f"{d} " + " ".join(map(str, s)) + f" {rng.randint(1, 50)}\n")
    return ("nmfsimtop", f"{T}/{name}.txt", [])
def cub(name, shape, periodic=False):
    A = np.random.default_rng(1).integers(1, 200, size=shape)
    from scipy.ndimage import gaussian_filter
    A = 1 + (gaussian_filter(A.astype(float), 3) * 10).astype(int)
    with open(f"{T}/{name}.txt", "w") as f:
        f.write("\n".join(map(str, [len(shape)] + [(-s if periodic else s) for s in shape])) + "\n"); np.savetxt(f, A.flatten(order="F"), fmt="%d")
    return ("cubtop", f"{T}/{name}.txt", [])
def pts(name, n, dim, maxdim):
    X = np.random.default_rng(2).normal(size=(n, dim)); np.savetxt(f"{T}/{name}.txt", X)
    return ("points", f"{T}/{name}.txt", ["--engine", "explicit", "--maxdim", str(maxdim), "--coeff", "2"])

cases = {
    "LM Y2 n=40 (all triangles)": lm("lm40", 40, 2, 1.0),
    "LM Y2 n=70 (55%)": lm("lm70", 70, 2, 0.55),
    "LM Y2 n=100 (40%)": lm("lm100", 100, 2, 0.4),
    "LM Y3 n=25 (all tetrahedra)": lm("lm3_25", 25, 3, 1.0),
    "random simplicial n=40, 3000 facets dim<=4": randsimp("rs", 40, 3000, 4),
    "cubical 1000^2": cub("c2", (1000, 1000)),
    "cubical 80^3": cub("c3", (80, 80, 80)),
    "cubical torus 512^2": cub("ct", (512, 512), True),
    "Rips explicit 120 pts R^4 H<=2": pts("r1", 120, 4, 2),
    "Rips explicit 60 pts R^8 H<=3": pts("r2", 60, 8, 3),
}
rows = []
for name, (fmt, path, extra) in cases.items():
    res = {}; outs = {}
    for d in ("cohomology", "homology", "forward"):
        pre = f"{T}/o_{d}"
        r = subprocess.run([P6, fmt, path, pre, "--stats", "--direction", d] + (extra or ["--coeff", "2"]), capture_output=True, text=True, timeout=900)
        m = re.search(r"(field run|Q run) ([0-9.]+) s", r.stderr); res[d] = float(m.group(2)) if m else None
        outs[d] = [sorted(open(f"{pre}_{k}.txt").read().splitlines()) if os.path.exists(f"{pre}_{k}.txt") else [] for k in range(6)]; [os.remove(f"{pre}_{k}.txt") for k in range(6) if os.path.exists(f"{pre}_{k}.txt")]
    same = outs["cohomology"] == outs["homology"] == outs["forward"]
    best = min(res, key=lambda k: res[k])
    rows.append(dict(case=name, seconds=res, agree=same, best=best)); print(json.dumps(rows[-1]), flush=True)
json.dump(rows, open(sys.argv[1] if len(sys.argv) > 1 else "direction.json", "w"), indent=1)
