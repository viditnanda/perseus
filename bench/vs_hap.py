"""Group homology and mod-p cohomology rings: perseus pcgroup against HAP's best routines (GAP).

Integral homology H_n(G; Z), n <= D: HAP is timed with GroupHomology(G, D) (HAP chooses the method)
and, for nilpotent groups, with ResolutionNilpotentGroup(G, D + 1) followed by Homology(TensorWithIntegers(R), k);
the faster of the two is reported. ResolutionNilpotentGroup builds a resolution of the same size as Perseus's
Anick resolution (C(n+k-1, k-1) generators in degree n for k cyclic factors; HAP may use fewer, larger cyclic
factors), so the comparison is between implementations of the same idea.
Mod-2 rings, degrees <= D: HAP's ModPCohomologyRing(G, D) and ModPRingGenerators (structure constants and minimal
generators, no relations), against perseus --ring 2 --no-homology (dims, all products, minimal generators).
The answers must agree (homology groups; degrees of the minimal generators); the script fails otherwise.

Times: HAP's Runtime() around the computation (GAP's start-up, package loading and SmallGroup excluded); Perseus's
own "seconds" (reading and validating the table excluded). Both are the best of three runs, except HAP runs over
REPEAT_BELOW = 10 s, which are run once. GroupHomology is timed as one call at the top degree, GroupHomology(G, D);
the lower degrees are computed afterwards, untimed, for the cross-check (each call builds its own resolution, and the
single call omits the lower-degree Smith forms, which favours HAP slightly). Each HAP process is limited to LIMIT
seconds of wall clock (default 300), start-up and cross-check included; a run is recorded as "timeout" only when that
limit stops it during the timed computation, and as "error", with the tail of GAP's output, when GAP fails. If the
limit falls in GroupHomology's untimed cross-check, the timed call's figure stands and only H_D is compared.

    python3 bench/vs_hap.py [--quick] [--out bench/vs_hap.json]

Needs: ./perseus (make), GAP with HAP on the PATH (run with GAP 4.12.1, HAP 1.62). Groups are SmallGroup(n, i),
written for perseus by tools/export_pcgroup.g. Not compared: the Green-King package for SageMath, the reference
for mod-p cohomology rings of p-groups.
"""
import json, os, re, subprocess, sys, tempfile, collections, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
PERSEUS = os.environ.get("PERSEUS", os.path.join(ROOT, "perseus"))
EXPORT = os.path.join(ROOT, "tools", "export_pcgroup.g")
TMP = tempfile.mkdtemp()
LIMIT = int(os.environ.get("LIMIT", "300"))
REPEAT_BELOW = 10.0

HOMOLOGY = [((16, 14), 6), ((16, 9), 8), ((24, 12), 6), ((27, 3), 6), ((32, 51), 6), ((32, 49), 6), ((64, 138), 6),
            ((64, 138), 8), ((128, 1000), 6), ((128, 2000), 6), ((243, 2), 6)]
RING = [((32, 49), 4), ((64, 138), 4), ((128, 1000), 4), ((64, 138), 8), ((128, 1000), 8)]
QUICK_H = [((16, 14), 6), ((27, 3), 6), ((64, 138), 6)]
QUICK_R = [((32, 49), 4)]


SIDE = os.path.join(TMP, "timed.txt")


def gap_once(body):
    """One GAP run. Returns {"status": "ok" | "timeout" | "error", "runtime_s": HAP's Runtime() around the timed part
    (ok only), "wall_s": wall clock of the whole process, "result": the printed result (ok only), "stderr": tail}.
    A body may write "MS <ms> RES <result> END" to SIDE (AppendTo closes the file, so it survives the process being
    killed) as soon as its timed part finishes; if the process then times out in its untimed cross-check, the timed
    figure is kept: status "timeout", stage "cross-check", with "runtime_s" and "result" from SIDE."""
    path = os.path.join(TMP, "s.g")
    open(path, "w").write('LoadPackage("hap");; SizeScreen([4096, 24]);;\n' + body + "\nQUIT;\n")
    if os.path.exists(SIDE): os.remove(SIDE)
    t0 = time.time()
    try:
        p = subprocess.run(["gap", "-q", "-o", "16g", path], capture_output=True, text=True, timeout=LIMIT)
    except subprocess.TimeoutExpired:
        r = {"status": "timeout", "wall_s": round(time.time() - t0, 1), "stage": "timed computation"}
        m = re.search(r"MS (\d+) RES (.*?) END", open(SIDE).read().replace("\\\n", "").replace("\n", " ")) if os.path.exists(SIDE) else None
        if m: r.update(stage="cross-check", runtime_s=int(m.group(1)) / 1000, result=json.loads(m.group(2).replace(" ", "")))
        return r
    wall = round(time.time() - t0, 3)
    out = p.stdout.replace("\\\n", "").replace("\n", " ")
    m = re.search(r"MS (\d+) RES (.*?) END", out)
    if not m:
        return {"status": "error", "wall_s": wall, "stderr": (p.stderr + p.stdout)[-400:]}
    return {"status": "ok", "runtime_s": int(m.group(1)) / 1000, "wall_s": wall, "result": json.loads(m.group(2).replace(" ", ""))}


def gap(body):
    """Best of three runs when the first takes under REPEAT_BELOW seconds of HAP's time; one run otherwise."""
    r = gap_once(body)
    if r["status"] == "ok" and r["runtime_s"] < REPEAT_BELOW:
        for _ in range(2):
            r2 = gap_once(body)
            if r2["status"] == "ok" and r2["runtime_s"] < r["runtime_s"]: r = r2
        r["runs"] = 3
    elif r["status"] == "ok":
        r["runs"] = 1
    return r


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
        runs = {}
        # GroupHomology: one call at the top degree is timed; the lower degrees are computed afterwards, untimed,
        # only for the cross-check (each GroupHomology call builds its own resolution)
        r = gap(G + ' t := Runtime();; top := GroupHomology(G, %d);; ms := Runtime() - t;; AppendTo("%s", "MS ", ms, " RES ", top, " END\\n");; H := Concatenation(List([1..%d], k -> GroupHomology(G, k)), [top]);; '
                'Print("MS ", ms, " RES ", H, " END\\n");' % (D, SIDE, D - 1))
        if r.get("stage") == "cross-check":
            # the timed call finished; only its top degree can be compared, and its time stands
            if prime_powers(r["result"]) != ours[-1]: raise SystemExit("MISMATCH (GroupHomology, top degree) on SmallGroup%s" % (g,))
            r = {"status": "ok", "runtime_s": r["runtime_s"], "wall_s": r["wall_s"], "runs": 1,
                 "note": "cross-check of degrees < %d stopped by the time limit; H_%d compared" % (D, D)}
        runs["GroupHomology"] = r
        if r["status"] == "ok" and "result" in r and [prime_powers(x) for x in r["result"]] != ours: raise SystemExit("MISMATCH (GroupHomology) on SmallGroup%s" % (g,))
        r = gap(G + ' if IsNilpotentGroup(G) then t := Runtime();; R := ResolutionNilpotentGroup(G, %d);; H := List([1..%d], k -> Homology(TensorWithIntegers(R), k));; '
                'Print("MS ", Runtime() - t, " RES ", H, " END\\n"); else Print("MS 0 RES \\"not_nilpotent\\" END\\n"); fi;' % (D + 1, D))
        if r["status"] == "ok" and r["result"] == "not_nilpotent": r = {"status": "not applicable"}
        runs["ResolutionNilpotentGroup"] = r
        if r["status"] == "ok" and [prime_powers(x) for x in r["result"]] != ours: raise SystemExit("MISMATCH (ResolutionNilpotentGroup) on SmallGroup%s" % (g,))
        ok = {k: v["runtime_s"] for k, v in runs.items() if v["status"] == "ok"}
        save({"group": "SmallGroup(%d, %d)" % g, "task": "H_n(G; Z), n <= %d" % D, "perseus_s": round(pt, 4),
              "hap": {k: {kk: vv for kk, vv in v.items() if kk != "result"} for k, v in runs.items()},
              "hap_best_s": min(ok.values()) if ok else None, "hap_best": min(ok, key=ok.get) if ok else None,
              "anick_chains_top": pj["chains"][-1]})
    for g, D in (QUICK_R if quick else RING):
        f = export(g); pj, pt = perseus(f, D, ring=2)
        r = gap('G := SmallGroup(%d, %d);; t := Runtime();; A := ModPCohomologyRing(G, %d);; gens := ModPRingGenerators(A);; '
                'Print("MS ", Runtime() - t, " RES ", List(gens, x -> A!.degree(x)), " END\\n");' % (g[0], g[1], D))
        ours = pj["cohomology"]["generators_by_degree"]
        if r["status"] == "ok":
            c = collections.Counter(d for d in r["result"] if d > 0)
            if [c.get(n, 0) for n in range(1, D + 1)] != ours[1:]: raise SystemExit("MISMATCH (ring generators) on SmallGroup%s" % (g,))
        save({"group": "SmallGroup(%d, %d)" % g, "task": "mod-2 ring, degrees <= %d" % D, "perseus_s": round(pt, 4),
              "hap": {"ModPCohomologyRing + ModPRingGenerators": {kk: vv for kk, vv in r.items() if kk != "result"}},
              "hap_best_s": r.get("runtime_s"), "hap_best": "ModPCohomologyRing + ModPRingGenerators" if r["status"] == "ok" else None,
              "generators_by_degree": ours[1:]})


if __name__ == "__main__":
    main()
