#!/usr/bin/env python3
"""Time/memory race: Perseus 6 (Q certified, and F_2) vs Ripser and RedZeD on lower-triangular float32 inputs."""
import json, os, resource, subprocess, sys, time
P6 = os.environ.get("PERSEUS", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "perseus"))
RIP = os.environ.get("RIPSER", "/home/claude/ext/ripser/ripser"); RZ = os.environ.get("REDZED", "/home/claude/build/redzed")
B = os.environ.get("BENCH", "/home/claude/bench/bin")
cases = [("sphere", 2, None), ("o3_1024", 3, 1.8), ("o3_4096", 1, 1.4), ("o3_4096", 3, 1.4), ("random16", 1, None),
         ("random16", 7, None), ("fractal", 1, None), ("fractal", 2, None), ("dragon", 1, None), ("randdist180", 1, None)]
def timed(cmd):
    t = time.time()
    r = subprocess.run(["python3", "-c", "import resource,subprocess,sys; subprocess.run(sys.argv[1:], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL); print(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)"] + cmd, capture_output=True, text=True)
    return time.time() - t, int(r.stdout.strip()) / 1024
rows = []
for f, dim, thr in cases:
    T = ["--threshold", str(thr)] if thr else []
    path = f"{B}/{f}.bin"
    q = timed([P6, "bin32", path, "--maxdim", str(dim), "--print", "none"] + T)
    f2 = timed([P6, "bin32", path, "--maxdim", str(dim), "--print", "none", "--coeff", "2"] + T)
    rp = timed([RIP, "--format", "binary", "--dim", str(dim)] + T + [path])
    rz = timed([RZ, path, str(dim), "--lower", "--barcode"] + T)
    row = dict(case=f"{f} H<={dim}", perseus_Q_s=round(q[0], 2), perseus_Q_MB=round(q[1]), perseus_F2_s=round(f2[0], 2), perseus_F2_MB=round(f2[1]),
               ripser_s=round(rp[0], 2), ripser_MB=round(rp[1]), redzed_s=round(rz[0], 2), redzed_MB=round(rz[1]))
    rows.append(row); print(row, flush=True)
json.dump(rows, open(sys.argv[1] if len(sys.argv) > 1 else "race.json", "w"), indent=1)
