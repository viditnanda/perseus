"""Group homology and mod-p cohomology rings: perseus pcgroup against HAP's best routines (GAP).

Integral homology H_n(G; Z), n <= D: HAP is timed with GroupHomology(G, k) for k = 1..D (HAP chooses the method)
and, for nilpotent groups, with ResolutionNilpotentGroup(G, D + 1) followed by Homology(TensorWithIntegers(R), k);
the faster of the two is reported. ResolutionNilpotentGroup builds a resolution of the same size as Perseus's
Anick resolution (C(n+k-1, k-1) generators in degree n for k cyclic factors; HAP may use fewer, larger cyclic
factors), so the comparison is between implementations of the same idea.
Mod-2 rings, degrees <= D: HAP's ModPCohomologyRing(G, D) and ModPRingGenerators (structure constants and minimal
generators, no relations), against perseus --ring 2 --no-homology (dims, all products, minimal generators).
The answers must agree (homology groups; degrees of the minimal generators); the script fails otherwise.

Times: HAP's Runtime() around the computation (GAP's start-up and package loading excluded); Perseus's own
"seconds" (reading the table excluded), best of three. Each HAP run is limited to LIMIT seconds (default 300).

    python3 bench/vs_hap.py [--quick] [--out bench/vs_hap.json]

Needs: ./perseus (make), GAP with HAP on the PATH (run with GAP 4.12.1, HAP 1.62). Groups are SmallGroup(n, i),
written for perseus by tools/export_pcgroup.g. Not compared: the Green-King package for SageMath, the reference
for mod-p cohomology rings of p-groups.
"""
import json, os, re, subprocess, sys, tempfile, collections

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
PERSEUS = os.environ.get("PERSEUS", os.path.join(ROOT, "perseus"))
EXPORT = os.path.join(ROOT, "tools", "export_pcgroup.g")
TMP = tempfile.mkdtemp()
LIMIT = int(os.environ.get("LIMIT", "300"))

HOMOLOGY = [((16, 14), 6), ((16, 9), 8), ((24, 12), 6), ((27, 3), 6), ((32, 51), 6), ((32, 49), 6), ((64, 138), 6),
            ((64, 138), 8), ((128, 1000), 6), ((128, 2000), 6), ((243, 2), 6)]
RING = [((32, 49), 4), ((64, 138), 4), ((128, 1000), 4), ((64, 138), 8), ((128, 1000), 8)]
QUICK_H = [((16, 14), 6), ((27, 3), 6), ((64, 138), 6)]
QUICK_R = [((32, 49), 4)]


def gap(body):
    path = os.path.join(TMP, "s.g")
    open(path, "w").write('LoadPackage("hap");; SizeScreen([4096, 24]);;\n' + body + "\nQUIT;\n")
    try:
        p = subprocess.run(["gap", "-q", "-o", "16g", path], capture_output=True, text=True, timeout=LIMIT + 60)
    except subprocess.TimeoutExpired:
        return None, None
    out = p.stdout.replace("\\\n", "").replace("\n", " ")
    m = re.search(r"MS (\d+) RES (.*?) END", out)
    if not m: return None, None
    t = int(m.group(1)) / 1000
    return (t if t <= LIMIT else None), json.loads(m.group(2).replace(" ", ""))


def export(g):
    f = os.path.join(TMP, "g%d_%d.txt" % g)
    if not os.path.exists(f):
        subprocess.run(["gap", "-q"], input='Read("%s");; WritePerseusPcGroup(SmallGroup(%d, %d), "%s");; QUIT;\n' % (EXPORT, g[0], g[1], f), text=True, capture_output=True)
        if not os.path.exists(f): raise RuntimeError("export failed for SmallGroup%s" % (g,))
    return f


def perseus(f, D, ring=0):
    best, j = None, None
    for _ in range(3):
        args = [PERSEUS, "pcgroup", f, "--maxdim", str(D), "--print", "json"] + (["--ring", str(ring), "--no-homology"] if ring else [])
        p = subprocess.run(args, capture_output=True, text=True)
        if p.returncode: raise RuntimeError(p.stderr)
        j = json.loads(p.stdout); best = j["seconds"] if best is None else min(best, j["seconds"])
    return j, best


def prime_powers(lst):
    out = []
    for x in lst:
        p = 2
        while x > 1:
            if x % p == 0:
                q = 1
                while x % p == 0: x //= p; q *= p
                out.append(q)
            p += 1
    return sorted(out)


def main():
    quick = "--quick" in sys.argv
    outpath = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else os.path.join(HERE, "vs_hap.json")
    rows = []
    def save(row):
        rows.append(row); print(json.dumps(row), flush=True)
        json.dump({"limit_s": LIMIT, "rows": rows}, open(outpath, "w"), indent=1)
    for g, D in (QUICK_H if quick else HOMOLOGY):
        f = export(g); pj, pt = perseus(f, D)
        ours = [prime_powers(h["torsion"]) for h in pj["homology"][1:]]
        G = "G := SmallGroup(%d, %d);;" % g
        times = {}
        t, H = gap(G + ' t := Runtime();; H := List([1..%d], k -> GroupHomology(G, k));; Print("MS ", Runtime() - t, " RES ", H, " END\\n");' % D)
        times["GroupHomology"] = t
        if t is not None and [prime_powers(x) for x in H] != ours: raise SystemExit("MISMATCH (GroupHomology) on SmallGroup%s" % (g,))
        t, H = gap(G + ' if IsNilpotentGroup(G) then t := Runtime();; R := ResolutionNilpotentGroup(G, %d);; H := List([1..%d], k -> Homology(TensorWithIntegers(R), k));; '
                   'Print("MS ", Runtime() - t, " RES ", H, " END\\n"); fi;' % (D + 1, D))
        times["ResolutionNilpotentGroup"] = t
        if t is not None and [prime_powers(x) for x in H] != ours: raise SystemExit("MISMATCH (ResolutionNilpotentGroup) on SmallGroup%s" % (g,))
        done = {k: v for k, v in times.items() if v is not None}
        save({"group": "SmallGroup(%d, %d)" % g, "task": "H_n(G; Z), n <= %d" % D, "perseus_s": round(pt, 4), "hap_times_s": times,
              "hap_best_s": min(done.values()) if done else None, "hap_best": min(done, key=done.get) if done else None,
              "anick_chains_top": pj["chains"][-1]})
    for g, D in (QUICK_R if quick else RING):
        f = export(g); pj, pt = perseus(f, D, ring=2)
        t, degs = gap('G := SmallGroup(%d, %d);; t := Runtime();; A := ModPCohomologyRing(G, %d);; gens := ModPRingGenerators(A);; '
                      'Print("MS ", Runtime() - t, " RES ", List(gens, x -> A!.degree(x)), " END\\n");' % (g[0], g[1], D))
        ours = pj["cohomology"]["generators_by_degree"]
        if t is not None:
            c = collections.Counter(d for d in degs if d > 0)
            if [c.get(n, 0) for n in range(1, D + 1)] != ours[1:]: raise SystemExit("MISMATCH (ring generators) on SmallGroup%s" % (g,))
        save({"group": "SmallGroup(%d, %d)" % g, "task": "mod-2 ring, degrees <= %d" % D, "perseus_s": round(pt, 4),
              "hap_best_s": t, "hap_best": "ModPCohomologyRing + ModPRingGenerators", "generators_by_degree": ours[1:]})


if __name__ == "__main__":
    main()
