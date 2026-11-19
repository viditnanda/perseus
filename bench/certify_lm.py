#!/usr/bin/env python3
"""Where the time goes when certifying a torsion-heavy complex: Linial-Meshulam processes near the threshold.

Instances: all edges of K_n at time 1, then the first m triangles of a uniformly random order at times
2, 3, ...; n = 70, m = 2400 and n = 90, m = 4000 (near c_2 * n^2 / 6 triangles, where H_1 dies and torsion
appears). Generated deterministically into /tmp/certify_lm. For each instance, the F_2 time, and the default
Q certificate with its phases from PERSEUS_DEBUG: the Q run, reconstruction and factoring of the pivot values,
the reruns for big candidates (primes >= 2^63, modulo multi-word N) and for small ones.

usage: bench/certify_lm.py [out.json] [--threads n]
"""
import json, os, re, subprocess, sys, time
from itertools import combinations
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
P6 = os.environ.get("PERSEUS", os.path.join(HERE, "..", "perseus"))
T = "/tmp/certify_lm"
CASES = [(70, 2400, 70), (90, 4000, 90)]      # n, triangles, seed


def make(n, m, seed):
    path = f"{T}/lm{n}_{m}.txt"
    if not os.path.exists(path):
        rng = np.random.default_rng(seed)
        tris = list(combinations(range(n), 3)); idx = rng.permutation(len(tris))[:m]
        with open(path, "w") as fh:
            for a, b in combinations(range(n), 2): fh.write(f"{a} {b} 1\n")
            for t, i in enumerate(idx): fh.write(" ".join(map(str, tris[i])) + f" {t + 2}\n")
    return path


def run(args, env=None):
    t = time.perf_counter()
    p = subprocess.run([P6] + args, capture_output=True, text=True, env=env)
    return time.perf_counter() - t, p.stdout, p.stderr


def main():
    out = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else os.path.join(HERE, "certify_lm.json")
    threads = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--threads=")), str(os.cpu_count()))
    os.makedirs(T, exist_ok=True)
    rows = []
    for n, m, seed in CASES:
        path = make(n, m, seed)
        row = {"case": f"Linial-Meshulam n={n}, first {m} triangles", "threads": int(threads)}
        row["F2_s"], _, _ = run(["simplices", path, "--coeff", "2", "--print", "none"])
        env = dict(os.environ, PERSEUS_DEBUG="1")
        row["Q_s"], so, se = run(["simplices", path, "--print", "json", "--threads", threads], env)
        j = json.loads(so)
        row["bad_primes"] = [str(b["p"]) for b in j["bad_primes"]]
        row["complete"] = j["complete"]
        g = lambda pat: re.search(pat, se)
        m1 = g(r"Q run ([0-9.]+) s, (\d+) non-unit pivots, largest (\d+) bits"); row["q_run_s"], row["nonunit"], row["largest_bits"] = float(m1[1]), int(m1[2]), int(m1[3])
        m2 = g(r"factoring ([0-9.]+) s, (\d+) \+ (\d+) candidates"); row["factoring_s"], row["small_candidates"], row["big_candidates"] = float(m2[1]), int(m2[2]), int(m2[3])
        m3 = g(r"big candidates ([0-9.]+) s, (\d+) runs"); row["big_s"], row["big_runs"] = float(m3[1]), int(m3[2])
        m4 = g(r"small candidates ([0-9.]+) s wall, ([0-9.]+) s summed"); row["small_wall_s"], row["small_cpu_s"] = float(m4[1]), float(m4[2])
        rows.append(row); print(json.dumps(row), flush=True)
    json.dump(rows, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
