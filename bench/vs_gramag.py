"""Magnitude homology: perseus (integral, certified) against gramag (ranks over F_2).

For every graph, both programs compute MH_{k,l} for l <= lmax (all pairs, or the pair (0^, 1^) of a
Kaneta-Yoshinaga graph). The F_2 ranks implied by perseus's integral groups (universal coefficients:
rank_F2 MH_k = free_k + #{even torsion in degree k} + #{even torsion in degree k-1}) must equal gramag's
ranks; the script fails loudly otherwise. Each run is a separate process, timed, with its peak RSS.

    python3 bench/vs_gramag.py [--quick] [--out bench/vs_gramag.json]

Needs: ./perseus (make), bench/mag_families (built by this script), pip install gramag. Each run is limited to
MEM_MB megabytes of address space (default 3500); a run that exceeds it is reported as failed.

What gramag computes (version 0.4.1): ranks of MH_{k,l} over Z/2 (its documentation: representatives are
"Z_2-linear combinations of paths"; the reduction is lophat's, which works over Z/2). It reports no torsion, so
the comparison is of F_2 ranks; the torsion column is Perseus's alone.
"""
import json, os, random, resource, subprocess, sys, tempfile, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
PERSEUS = os.environ.get("PERSEUS", os.path.join(ROOT, "perseus"))
FAM = os.path.join(tempfile.gettempdir(), "perseus_mag_families")
TMP = tempfile.mkdtemp()
THREADS = os.environ.get("THREADS", "2")


MEM_MB = int(os.environ.get("MEM_MB", "3500"))     # address-space limit for each run (the machine is shared)


class Failed(Exception):
    pass


def measured(cmd):
    """run cmd in a fresh wrapper process under the memory limit; return (stdout, seconds, peak MB of cmd)"""
    wrap = ("import resource, subprocess, sys, time; t = time.time(); "
            f"lim = lambda: resource.setrlimit(resource.RLIMIT_AS, ({MEM_MB} << 20, {MEM_MB} << 20)); "
            "r = subprocess.run(sys.argv[1:], capture_output=True, text=True, preexec_fn=lim); dt = time.time() - t; "
            "sys.stderr.write(r.stderr); "
            "print(r.returncode, dt, resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss); print(r.stdout, end='')")
    out = subprocess.run([sys.executable, "-c", wrap] + cmd, capture_output=True, text=True)
    head, _, body = out.stdout.partition("\n")
    rc, dt, rss = head.split()
    if int(rc) != 0: raise Failed(f"exit {rc} after {float(dt):.1f}s (memory limit {MEM_MB} MB): {out.stderr.strip()[-200:]}")
    return body, float(dt), int(rss) / 1024.0


GRAMAG = r'''
import json, sys, gramag
edges, lmax, pairs = json.load(open(sys.argv[1]))
g = gramag.MagGraph([tuple(e) for e in edges])
g.populate_paths(l_max=lmax)
r = g.rank_homology(node_pairs=[tuple(p) for p in pairs]) if pairs else g.rank_homology()
print(json.dumps(r))
'''


def gramag_ranks(edges, directed, lmax, pairs):
    E = edges if directed else edges + [(b, a) for a, b in edges]
    path = os.path.join(TMP, "g.json"); json.dump([E, lmax, pairs], open(path, "w"))
    script = os.path.join(TMP, "gramag_run.py"); open(script, "w").write(GRAMAG)
    out, dt, mb = measured([sys.executable, script, path])
    r = json.loads(out)
    ranks = {}
    for l, row in enumerate(r):
        for k, x in enumerate(row):
            if x: ranks[(k, l)] = x
    return ranks, dt, mb


def perseus_groups(edges, directed, lmax, pairs):
    path = os.path.join(TMP, "g.txt")
    with open(path, "w") as f:
        for a, b in edges: f.write(f"{a} {b}\n")
    cmd = [PERSEUS, "magnitude", path, "--lmax", str(lmax), "--print", "json", "--threads", THREADS]
    if directed: cmd.append("--directed")
    if pairs: cmd += ["--pairs", ",".join(f"{s},{t}" for s, t in pairs)]
    out, dt, mb = measured(cmd)
    j = json.loads(out)
    return {(g["k"], g["l"]): g for g in j["groups"]}, dt, mb


def f2_ranks(G):
    ev = lambda g: sum(m for q, m in g["torsion"].items() if int(q) % 2 == 0) if g else 0
    r = {}
    for (k, l), g in G.items():
        r[(k, l)] = g["free"] + ev(g) + ev(G.get((k - 1, l)))
    return {key: v for key, v in r.items() if v}


def ky(case):
    path = os.path.join(TMP, case + ".edges")
    subprocess.run([FAM, case, "0", "pair", path], check=True, capture_output=True)
    lines = open(path).read().split("\n")
    n = int(lines[0]); E = [tuple(map(int, x.split())) for x in lines[1:] if x.strip()]
    return n, E


def main():
    quick = "--quick" in sys.argv
    outp = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else os.path.join(HERE, "vs_gramag.json")
    subprocess.run(["g++", "-std=c++17", "-O2", "-pthread", "-I" + os.path.join(ROOT, "include"),
                    os.path.join(HERE, "mag_families.cpp"), "-o", FAM], check=True)
    rng = random.Random(7)
    cases = []
    cyc = lambda n: [(i, (i + 1) % n) for i in range(n)]
    cases.append(("C_7", cyc(7), False, 7, []))
    cases.append(("C_12", cyc(12), False, 8 if quick else 10, []))
    cases.append(("directed C_6", cyc(6), True, 8, []))
    grid = [(r * 4 + c, r * 4 + c + 1) for r in range(4) for c in range(3)] + [(r * 4 + c, r * 4 + c + 4) for r in range(3) for c in range(4)]
    cases.append(("grid 4x4", grid, False, 5 if quick else 6, []))
    er = sorted({(i, j) for i in range(20) for j in range(i + 1, 20) if rng.random() < 0.2})
    cases.append(("G(20, 0.2)", er, False, 4 if quick else 5, []))
    dr = sorted({(i, j) for i in range(16) for j in range(16) if i != j and rng.random() < 0.18})
    cases.append(("directed G(16, 0.18)", dr, True, 4 if quick else 5, []))
    for case, lmax, allpairs in (("rp2", 6, True), ("rp2_6", 5, True), ("rp2", 8, False), ("rp3", 7, False), ("lens3_2", 7, False)):
        n, E = ky(case)
        pairs = [] if allpairs else [(n - 2, n - 1)]
        if quick and lmax > 6: lmax = 6
        cases.append((f"KY {case} ({'all pairs' if allpairs else 'pair (0,1)'})", E, False, lmax, pairs))
    rows, bad = [], 0
    for name, E, directed, lmax, pairs in cases:
        G, tp, mp = perseus_groups(E, directed, lmax, pairs)
        try:
            ranks, tg, mg = gramag_ranks(E, directed, lmax, pairs)
        except Failed as e:
            ranks, tg, mg = None, None, None
            print(f"    gramag failed: {e}", flush=True)
        mine = f2_ranks(G)
        ok = mine == ranks if ranks is not None else None
        bad += ok is False
        tors = sorted((k, l, {int(q): m for q, m in g["torsion"].items()}) for (k, l), g in G.items() if g["torsion"])
        gens = sum(g["generators"] for g in G.values())
        n = 1 + max(max(e) for e in E)
        row = dict(graph=name, vertices=n, edges=len(E), lmax=lmax, generators=gens, agree_over_F2=ok,
                   perseus_s=round(tp, 3), perseus_mb=round(mp, 1), gramag_s=None if tg is None else round(tg, 3),
                   gramag_mb=None if mg is None else round(mg, 1), torsion=tors)
        rows.append(row)
        gm = "  failed" + " " * 10 if tg is None else f"{tg:7.2f}s {mg:7.1f}MB"
        verdict = "n/a" if ok is None else ("agree" if ok else "DISAGREE")
        print(f"{name:32s} n={n:3d} lmax={lmax:2d} gens={gens:9d}  perseus {tp:7.2f}s {mp:7.1f}MB   gramag {gm}   "
              f"F2 ranks {verdict}  torsion {tors}", flush=True)
        if ok is False:
            for key in sorted(set(mine) | set(ranks)):
                if mine.get(key) != ranks.get(key): print("    ", key, "perseus", mine.get(key), "gramag", ranks.get(key))
    json.dump(rows, open(outp, "w"), indent=1)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
