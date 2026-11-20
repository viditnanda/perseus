#!/usr/bin/env python3
"""Where the time goes when certifying a torsion-heavy complex: Linial-Meshulam processes near the threshold.

Instances: all edges of K_n at time 1, then the first m triangles of a uniformly random order at times
2, 3, ...; n = 70, m = 2400 and n = 90, m = 4000 (near c_2 * n^2 / 6 triangles, where H_1 dies and torsion
appears). Generated deterministically into /tmp/certify_lm. For each instance, the F_2 time; the default Q
certificate (proved: exact rationals) with its phases from PERSEUS_DEBUG: the Q run, factoring of the pivot
values, the reruns for big candidates (primes >= 2^63, modulo multi-word N) and for small ones; and, for
comparison, --mode fast (random primes, not a proof), checked to find the same bad primes. Small candidates are
classified by batched Howell sweeps (the default for explicit complexes since 19 Nov 2026); the reruns they
replace are timed too (--classify-by rerun) and must give the same bad primes. Last, --torsion: the torsion of
H_1(K_t; Z) at every t, with its cost. A third instance has distinct, interleaved values (edges at random times,
each triangle after its edges), where almost every candidate is bad.

usage: bench/certify_lm.py [out.json] [--threads n]
"""
import json, os, re, subprocess, sys, time
from itertools import combinations
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
P6 = os.environ.get("PERSEUS", os.path.join(HERE, "..", "perseus"))
T = "/tmp/certify_lm"
CASES = [(70, 2400, 70, "tied"), (90, 4000, 90, "tied"), (70, 2400, 70, "distinct")]   # n, triangles, seed, values


def make(n, m, seed, kind="tied"):
    if kind == "distinct":
        path = f"{T}/lmd{n}_{m}.txt"
        if not os.path.exists(path):
            rng = np.random.default_rng(seed)
            E = {e: rng.random() for e in combinations(range(n), 2)}
            tris = list(combinations(range(n), 3)); idx = rng.permutation(len(tris))[:m]
            with open(path, "w") as fh:
                for (a, b), w in E.items(): fh.write(f"{a} {b} {w:.12f}\n")
                for i in idx:
                    a, b, c = tris[i]; me = max(E[(a, b)], E[(a, c)], E[(b, c)])
                    fh.write(f"{a} {b} {c} {max(rng.random() * 1.3, me + 1e-7 * (1 + rng.random())):.12f}\n")
        return path
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
    for n, m, seed, kind in CASES:
        path = make(n, m, seed, kind)
        row = {"case": f"Linial-Meshulam n={n}, first {m} triangles" + (", distinct interleaved values" if kind == "distinct" else ""), "threads": int(threads)}
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
        m4 = g(r"sweep classification used: (\d+) sweeps for (\d+) candidates in ([0-9.]+) s"); row["sweeps"], row["small_sweep_s"] = int(m4[1]), float(m4[3])
        row["rerun_Q_s"], so, se = run(["simplices", path, "--print", "json", "--threads", threads, "--classify-by", "rerun"], env)
        m5 = re.search(r"small candidates ([0-9.]+) s wall, ([0-9.]+) s summed over (\d+) reruns", se)
        row["small_rerun_wall_s"], row["small_rerun_cpu_s"], row["reruns"] = float(m5[1]), float(m5[2]), int(m5[3])
        row["rerun_same_bad_primes"] = json.loads(so)["bad_primes"] == j["bad_primes"]
        tt, so, _ = run(["simplices", path, "--print", "json", "--threads", threads, "--torsion"])
        jt = json.loads(so)["torsion"]
        row["torsion_s"], row["torsion_sweeps"] = jt["seconds"], jt["sweeps"]
        row["torsion_steps"] = [[s["value"], s["invariant_factors"]] for d in jt["degrees"] for s in d["steps"]][:12]
        row["mode"] = j["mode"]; row["proof"] = j.get("proof")
        row["fast_s"], so, _ = run(["simplices", path, "--print", "json", "--threads", threads, "--mode", "fast"])
        jf = json.loads(so); row["fast_bad_primes_agree"] = [str(b["p"]) for b in jf["bad_primes"]] == row["bad_primes"]
        rows.append(row); print(json.dumps(row), flush=True)
    json.dump(rows, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
