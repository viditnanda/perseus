"""Group homology and mod-p cohomology rings: perseus pcgroup against HAP's best routines (GAP).

Integral homology H_n(G; Z), n <= D: HAP is timed with GroupHomology(G, D) (HAP chooses the method)
and, for nilpotent groups, with ResolutionNilpotentGroup(G, D + 1) followed by Homology(TensorWithIntegers(R), k);
the faster of the two is reported. ResolutionNilpotentGroup builds a resolution of the same size as Perseus's
Anick resolution (C(n+k-1, k-1) generators in degree n for k cyclic factors; HAP may use fewer, larger cyclic
factors), so the comparison is between implementations of the same idea.
Mod-2 rings, degrees <= D: HAP's ModPCohomologyRing(G, D) and ModPRingGenerators (structure constants and minimal
generators, no relations), against perseus --ring 2 --no-homology --all-products (dims, all products of basis classes,
minimal generators), which is the same output. Perseus's default, the products of the generators with every class
(enough to find the generators), is timed too, as perseus_generators_only_s.
The answers must agree (homology groups; degrees of the minimal generators); the script fails otherwise.

Times: HAP's Runtime() around the computation (GAP's start-up, package loading and SmallGroup excluded); Perseus's
own "seconds" (reading and validating the table excluded). Both are the best of three runs, except HAP runs over
REPEAT_BELOW = 10 s, which are run once. GroupHomology is timed as one call at the top degree, GroupHomology(G, D);
the lower degrees are computed afterwards, untimed, for the cross-check (each call builds its own resolution, and the
single call omits the lower-degree Smith forms, which favours HAP slightly). Each HAP process is limited to LIMIT
seconds of wall clock (default 300), start-up and cross-check included; a run is recorded as "timeout" only when that
limit stops it during the timed computation, and as "error", with the tail of GAP's output, when GAP fails. If the
limit falls in GroupHomology's untimed cross-check, the timed call's figure stands and only H_D is compared.

    python3 bench/vs_hap.py [--quick] [--rings | --sylow | --series] [--out bench/vs_hap.json]   (rerun only those rows)

Series rows: persistence of mod-p homology along a quotient series (perseus pcseries, export by
tools/export_pcseries.g timed and added as for the Sylow rows) against HAP's PersistentHomologyOfQuotientGroupSeries
(S, n), which gives degree n only; Perseus computes every degree <= n. HAP 1.62 needs a one-line fix for this routine
(see hap_patch), applied in the GAP session only.

Sylow rows: any finite group, perseus sylow against HAP, which starts from G itself. Perseus starts from the Sylow
data that tools/export_sylow.g computes in GAP (Sylow subgroups, double cosets, intersections), so that export is
timed too (export_s, GAP's Runtime(), best of three) and perseus_total_s = perseus_s + export_s is the figure to
compare. Homology against GroupHomology(G, D), timed and cross-checked as above; mod-p rings of groups that are not
p-groups against ModPCohomologyRing(G, p, D) and ModPRingGenerators, with perseus sylow --all-products (the same
output) and generator degrees cross-checked.

Needs: ./perseus (make), GAP with HAP on the PATH (run with GAP 4.12.1, HAP 1.62). Groups are SmallGroup(n, i),
written for perseus by tools/export_pcgroup.g, which chooses a short pcgs in GAP; that export is timed like the Sylow
and series exports (export_s, best of three) and perseus_total_s = perseus_s + export_s. Not compared: the Green-King package for SageMath, the reference
for mod-p cohomology rings of p-groups.
"""
import collections, json, os, re, subprocess, sys, tempfile, time

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
# any finite group through its Sylow subgroups (perseus sylow): (GAP expression, name, D)
SYLOW = [("SmallGroup(48, 29)", "GL(2,3)", 6), ("SmallGroup(72, 40)", "(S3 x S3) : C2", 6), ("SmallGroup(96, 64)", "(C4 x C4) : S3", 6),
         ("SmallGroup(96, 64)", "(C4 x C4) : S3", 10), ("AlternatingGroup(5)", "A5", 8), ("PSL(2, 7)", "PSL(2,7)", 8),
         ("AlternatingGroup(6)", "A6", 8), ("AlternatingGroup(7)", "A7", 8), ("MathieuGroup(11)", "M11", 8),
         ("MathieuGroup(12)", "M12", 6), ("PSL(3, 4)", "PSL(3,4)", 6)]
# mod-p rings of groups that are not p-groups, through Sylow subgroups: (GAP expression, name, p, D)
SYLOW_RING = [("AlternatingGroup(5)", "A5", 2, 8), ("PSL(2, 7)", "PSL(2,7)", 2, 8), ("SymmetricGroup(4)", "S4", 3, 8),
              ("AlternatingGroup(6)", "A6", 3, 6), ("MathieuGroup(11)", "M11", 2, 6), ("MathieuGroup(11)", "M11", 2, 8),
              ("MathieuGroup(12)", "M12", 2, 6)]
# persistence along a series (perseus pcseries): (GAP group, series, name, n)
SERIES = [("SmallGroup(64, 138)", "LowerCentralSeries(G)", "SmallGroup(64,138), lower central", 6),
          ("SmallGroup(64, 138)", "LowerCentralSeries(G)", "SmallGroup(64,138), lower central", 8),
          ("DihedralGroup(128)", "LowerCentralSeries(G)", "D128, lower central", 8),
          ("SmallGroup(128, 1000)", "PCentralSeries(G, 2)", "SmallGroup(128,1000), lower 2-central", 6),
          ("SmallGroup(243, 2)", "LowerCentralSeries(G)", "SmallGroup(243,2), lower central", 6),
          ("SmallGroup(243, 26)", "LowerCentralSeries(G)", "SmallGroup(243,26), lower central", 6)]
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


EXPORT_S = {}


def export(g):
    """Write SmallGroup(g) with tools/export_pcgroup.g (a short pcgs, and GAP's as the alternative). The export is timed
    like the Sylow and series exports (GAP's Runtime() around WritePerseusPcGroup, best of three), since choosing the
    short pcgs is GAP-side work; EXPORT_S[g] holds the time."""
    f = os.path.join(TMP, "g%d_%d.txt" % g)
    if g not in EXPORT_S:
        best = None
        for _ in range(3):
            p = subprocess.run(["gap", "-q"], input='Read("%s");; G := SmallGroup(%d, %d);; t := Runtime();; WritePerseusPcGroup(G, "%s");; Print("MS ", Runtime() - t, "\\n");; QUIT;\n'
                               % (EXPORT, g[0], g[1], f), text=True, capture_output=True)
            m = re.search(r"MS (\d+)", p.stdout)
            if not m or not os.path.exists(f): raise RuntimeError("export failed for SmallGroup%s: %s" % (g, p.stderr[-300:]))
            t = int(m.group(1)) / 1000; best = t if best is None else min(best, t)
        EXPORT_S[g] = best
    return f



def export_sylow(expr):
    """Write the Sylow data of G; returns (path, seconds). The time is GAP's Runtime() around WritePerseusSylow
    (Sylow subgroups, double cosets, intersections, pcgs and tables, and writing the file; start-up excluded, as for
    HAP), best of three when under 10 s."""
    f = os.path.join(TMP, "sy_%s.txt" % re.sub(r"[^A-Za-z0-9]", "_", expr))
    if expr not in EXPORT_S:
        best = None
        for _ in range(3):
            p = subprocess.run(["gap", "-q"], input='Read("%s");; G := %s;; t := Runtime();; WritePerseusSylow(G, "%s");; Print("MS ", Runtime() - t, "\\n");; QUIT;\n'
                               % (os.path.join(ROOT, "tools", "export_sylow.g"), expr, f), text=True, capture_output=True)
            m = re.search(r"MS (\d+)", p.stdout)
            if not m or not os.path.exists(f): raise RuntimeError("Sylow export failed for " + expr + ": " + p.stderr[-500:])
            t = int(m.group(1)) / 1000; best = t if best is None else min(best, t)
            if best >= 10: break
        EXPORT_S[expr] = best
    return f, EXPORT_S[expr]


def hap_patch():
    """HAP 1.62's HomologyVectorSpace has the branches of its first test swapped (line 15 of
    lib/Homology/modularHomologyVectSpace.gi: C:=X for a chain map, Source(X) otherwise), which breaks
    PersistentHomologyOfQuotientGroupSeries. Returns GAP code loading a corrected copy, or "" if not needed."""
    src = "/usr/share/gap/pkg/HAP/lib/Homology/modularHomologyVectSpace.gi"
    bad = "if IsHapChainMap(X) then C:=X; else C:=Source(X); fi;"
    if not os.path.exists(src) or bad not in open(src).read(): return ""
    fixed = os.path.join(TMP, "hvs_fixed.gi")
    open(fixed, "w").write(open(src).read().replace(bad, "if IsHapChainMap(X) then C:=Source(X); else C:=X; fi;"))
    return ('MakeReadWriteGlobal("HomologyVectorSpace");; UnbindGlobal("HomologyVectorSpace");; '
            'DeclareGlobalFunction("HomologyVectorSpace");; Read("%s");; ' % fixed)


def export_series(gexpr, sexpr):
    key = gexpr + "|" + sexpr
    f = os.path.join(TMP, "ps_%s.txt" % re.sub(r"[^A-Za-z0-9]", "_", key))
    if key not in EXPORT_S:
        best = None
        for _ in range(3):
            p = subprocess.run(["gap", "-q"], input='Read("%s");; G := %s;; S := %s;; t := Runtime();; WritePerseusQuotientSeries(G, S, "%s");; Print("MS ", Runtime() - t, "\\n");; QUIT;\n'
                               % (os.path.join(ROOT, "tools", "export_pcseries.g"), gexpr, sexpr, f), text=True, capture_output=True)
            m = re.search(r"MS (\d+)", p.stdout)
            if not m or not os.path.exists(f): raise RuntimeError("series export failed for " + key + ": " + p.stderr[-500:])
            t = int(m.group(1)) / 1000; best = t if best is None else min(best, t)
            if best >= 10: break
        EXPORT_S[key] = best
    return f, EXPORT_S[key]


def perseus_sylow(f, D, extra=()):
    best = None
    for _ in range(3):
        p = subprocess.run([PERSEUS, "sylow", f, "--maxdim", str(D), "--print", "json"] + list(extra), capture_output=True, text=True)
        if p.returncode: raise RuntimeError(p.stderr)
        j = json.loads(p.stdout); best = j["seconds"] if best is None else min(best, j["seconds"])
    return j, best


def perseus(f, D, ring=0, extra=()):
    best, j = None, None
    for _ in range(3):
        args = [PERSEUS, "pcgroup", f, "--maxdim", str(D), "--print", "json"] + (["--ring", str(ring), "--no-homology"] if ring else []) + list(extra)
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
    rings_only = "--rings" in sys.argv
    sylow_only = "--sylow" in sys.argv
    series_only = "--series" in sys.argv
    if rings_only:   # rerun the ring rows only, keeping the other rows already in outpath
        rows = [r for r in json.load(open(outpath))["rows"] if not r["task"].startswith("mod-")]
    if sylow_only:   # rerun the Sylow rows only
        rows = [r for r in json.load(open(outpath))["rows"] if "Sylow" not in r["task"]]
    if series_only:  # rerun the series rows only
        rows = [r for r in json.load(open(outpath))["rows"] if "persistence" not in r["task"]]
    def save(row):
        rows.append(row); print(json.dumps(row), flush=True)
        json.dump({"limit_s": LIMIT, "rows": rows}, open(outpath, "w"), indent=1)
    patch = hap_patch()
    for gexpr, sexpr, name, n in (SERIES if series_only or not (quick or rings_only or sylow_only) else []):
        f, ex = export_series(gexpr, sexpr)
        best = None
        for _ in range(3):
            p = subprocess.run([PERSEUS, "pcseries", f, "--maxdim", str(n), "--print", "json"], capture_output=True, text=True)
            if p.returncode: raise RuntimeError(p.stderr)
            pj = json.loads(p.stdout); best = pj["seconds"] if best is None else min(best, pj["seconds"])
        ours = pj["degrees"][n]["rank"]
        r = gap(patch + 'G := %s;; S := %s;; t := Runtime();; P := PersistentHomologyOfQuotientGroupSeries(S, %d);; Print("MS ", Runtime() - t, " RES ", P, " END\\n");' % (gexpr, sexpr, n))
        if r["status"] == "ok":
            P = r["result"]; m = len(ours)
            if [[P[i][j] if j >= i else 0 for j in range(m)] for i in range(m)] != [[ours[i][j] if j >= i else 0 for j in range(m)] for i in range(m)]:
                raise SystemExit("MISMATCH (persistence) on " + name)
        save({"group": name, "task": "persistence of H_%d(-; F_p) along the series (Perseus: all degrees <= %d)" % (n, n), "perseus_s": round(best, 4),
              "export_s": round(ex, 4), "perseus_total_s": round(best + ex, 4), "quotients": len(ours),
              "hap": {"PersistentHomologyOfQuotientGroupSeries": {k: v for k, v in r.items() if k != "result"}},
              "hap_best_s": r.get("runtime_s") if r["status"] == "ok" else None, "hap_best": "PersistentHomologyOfQuotientGroupSeries" if r["status"] == "ok" else None,
              "rank": ours})
    for expr, name, D in (SYLOW if sylow_only or not (quick or rings_only or series_only) else []):
        f, ex = export_sylow(expr); pj, pt = perseus_sylow(f, D)
        ours = [prime_powers(h["torsion"]) for h in pj["homology"][1:]]
        r = gap('G := %s;; t := Runtime();; top := GroupHomology(G, %d);; ms := Runtime() - t;; AppendTo("%s", "MS ", ms, " RES ", top, " END\\n");; H := Concatenation(List([1..%d], k -> GroupHomology(G, k)), [top]);; '
                'Print("MS ", ms, " RES ", H, " END\\n");' % (expr, D, SIDE, D - 1))
        if r.get("stage") == "cross-check":
            if prime_powers(r["result"]) != ours[-1]: raise SystemExit("MISMATCH (GroupHomology, top degree) on " + name)
            r = {"status": "ok", "runtime_s": r["runtime_s"], "wall_s": r["wall_s"], "runs": 1, "note": "cross-check of lower degrees stopped by the time limit; H_%d compared" % D}
        elif r["status"] == "ok" and [prime_powers(x) for x in r["result"]] != ours: raise SystemExit("MISMATCH (GroupHomology) on " + name)
        save({"group": name, "task": "H_n(G; Z), n <= %d, via Sylow subgroups" % D, "perseus_s": round(pt, 4), "export_s": round(ex, 4), "perseus_total_s": round(pt + ex, 4),
              "hap": {"GroupHomology": {k: v for k, v in r.items() if k != "result"}},
              "hap_best_s": r.get("runtime_s") if r["status"] == "ok" else None, "hap_best": "GroupHomology" if r["status"] == "ok" else None,
              "answer": [",".join(map(str, x)) for x in ours]})
    for expr, name, p, D in (SYLOW_RING if sylow_only or not (quick or rings_only or series_only) else []):
        f, ex = export_sylow(expr)
        pj, pt = perseus_sylow(f, D, ["--ring", str(p), "--no-homology", "--all-products"])
        pg, ptg = perseus_sylow(f, D, ["--ring", str(p), "--no-homology"])
        ours = pj["cohomology"]["generators_by_degree"][1:]
        if pg["cohomology"]["generators_by_degree"][1:] != ours: raise SystemExit("MISMATCH (product modes) on " + name)
        r = gap('G := %s;; t := Runtime();; A := ModPCohomologyRing(G, %d, %d);; gens := ModPRingGenerators(A);; '
                'Print("MS ", Runtime() - t, " RES ", List(gens, x -> A!.degree(x)), " END\\n");' % (expr, p, D))
        if r["status"] == "ok":
            c = collections.Counter(d for d in r["result"] if d > 0)
            if [c.get(n, 0) for n in range(1, D + 1)] != ours: raise SystemExit("MISMATCH (ring generators) on " + name)
        save({"group": name, "task": "mod-%d ring, degrees <= %d, via Sylow subgroups" % (p, D), "perseus_s": round(pt, 4), "perseus_generators_only_s": round(ptg, 4),
              "export_s": round(ex, 4), "perseus_total_s": round(pt + ex, 4),
              "hap": {"ModPCohomologyRing": {k: v for k, v in r.items() if k != "result"}},
              "hap_best_s": r.get("runtime_s") if r["status"] == "ok" else None, "hap_best": "ModPCohomologyRing" if r["status"] == "ok" else None,
              "generators_by_degree": ours})
    for g, D in ([] if rings_only or sylow_only or series_only else QUICK_H if quick else HOMOLOGY):
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
        save({"group": "SmallGroup(%d, %d)" % g, "task": "H_n(G; Z), n <= %d" % D, "perseus_s": round(pt, 4), "export_s": round(EXPORT_S[g], 4), "perseus_total_s": round(pt + EXPORT_S[g], 4),
              "hap": {k: {kk: vv for kk, vv in v.items() if kk != "result"} for k, v in runs.items()},
              "hap_best_s": min(ok.values()) if ok else None, "hap_best": min(ok, key=ok.get) if ok else None,
              "anick_chains_top": pj["chains"][-1]})
    for g, D in ([] if sylow_only or series_only else QUICK_R if quick else RING):
        f = export(g); pj, pt = perseus(f, D, ring=2, extra=["--all-products"])
        pg, ptg = perseus(f, D, ring=2)
        if pg["cohomology"]["generators_by_degree"] != pj["cohomology"]["generators_by_degree"]: raise SystemExit("MISMATCH (product modes) on SmallGroup%s" % (g,))
        r = gap('G := SmallGroup(%d, %d);; t := Runtime();; A := ModPCohomologyRing(G, %d);; gens := ModPRingGenerators(A);; '
                'Print("MS ", Runtime() - t, " RES ", List(gens, x -> A!.degree(x)), " END\\n");' % (g[0], g[1], D))
        ours = pj["cohomology"]["generators_by_degree"]
        if r["status"] == "ok":
            c = collections.Counter(d for d in r["result"] if d > 0)
            if [c.get(n, 0) for n in range(1, D + 1)] != ours[1:]: raise SystemExit("MISMATCH (ring generators) on SmallGroup%s" % (g,))
        save({"group": "SmallGroup(%d, %d)" % g, "task": "mod-2 ring, degrees <= %d" % D, "perseus_s": round(pt, 4), "perseus_generators_only_s": round(ptg, 4), "export_s": round(EXPORT_S[g], 4), "perseus_total_s": round(pt + EXPORT_S[g], 4),
              "hap": {"ModPCohomologyRing + ModPRingGenerators": {kk: vv for kk, vv in r.items() if kk != "result"}},
              "hap_best_s": r.get("runtime_s"), "hap_best": "ModPCohomologyRing + ModPRingGenerators" if r["status"] == "ok" else None,
              "generators_by_degree": ours[1:]})


if __name__ == "__main__":
    main()
