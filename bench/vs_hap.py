"""Group homology and mod-p cohomology rings: perseus pcgroup against HAP (GAP).

For each group, both programs compute H_n(G; Z) for n <= D; the answers must agree (the script fails loudly
otherwise). For the mod-2 ring, perseus computes dims, all products and the minimal generators through degree D;
HAP's Mod2CohomologyRingPresentation(G, D) computes a presentation (generators and relations) through degree D.
Each run is a separate process, timed, with a wall-clock limit (LIMIT seconds, default 600).

    python3 bench/vs_hap.py [--quick] [--out bench/vs_hap.json]

Needs: ./perseus (make), GAP with the HAP package (`gap` on the PATH; this was run with GAP 4.12.1, HAP 1.62).
The groups are SmallGroup(n, i), written for perseus by tools/export_pcgroup.g.

What HAP does: ResolutionFiniteGroup builds a free ZG-resolution degree by degree by linear algebra over ZG,
with its rank kept small by a contracting homotopy (Ellis, "Computing group resolutions", JSC 2004). Perseus
uses the Anick resolution of the pcgs collection rewriting system (C(n+k-1, k-1) generators in degree n) and
cancels units over Z/q^(v+1) for each prime q | |G|. The mod-p ring is a different task in the two programs:
HAP also finds relations, which perseus does not; perseus's products are all products of basis elements.
The Green-King package for SageMath (the reference for mod-p rings of p-groups) is not compared here.
"""
import json, os, re, subprocess, sys, tempfile, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
PERSEUS = os.environ.get("PERSEUS", os.path.join(ROOT, "perseus"))
EXPORT = os.path.join(ROOT, "tools", "export_pcgroup.g")
TMP = tempfile.mkdtemp()
LIMIT = int(os.environ.get("LIMIT", "600"))

HOMOLOGY = [((16, 14), 6), ((16, 9), 8), ((24, 12), 6), ((27, 3), 6), ((32, 51), 6), ((32, 49), 6), ((64, 138), 6),
            ((64, 138), 8), ((128, 1000), 6), ((128, 2000), 6), ((243, 2), 6)]
RING = [((32, 49), 4), ((64, 138), 4), ((128, 1000), 4)]
QUICK_H = [((16, 14), 6), ((32, 51), 6), ((64, 138), 6)]
QUICK_R = [((32, 49), 4)]


def gap(script):
    path = os.path.join(TMP, "s.g"); open(path, "w").write(script + "\nQUIT;\n")
    t = time.time()
    try:
        p = subprocess.run(["gap", "-q", "-o", "16g", path], capture_output=True, text=True, timeout=LIMIT)
    except subprocess.TimeoutExpired:
        return None, LIMIT
    return p.stdout, time.time() - t


def export(g):
    f = os.path.join(TMP, "g%d_%d.txt" % g)
    if not os.path.exists(f):
        out, _ = gap('Read("%s");; WritePerseusPcGroup(SmallGroup(%d, %d), "%s");' % (EXPORT, g[0], g[1], f))
        if not os.path.exists(f): raise RuntimeError("export failed: %s" % out)
    return f


def perseus(f, D, ring=0):
    args = [PERSEUS, "pcgroup", f, "--maxdim", str(D), "--print", "json"] + (["--ring", str(ring), "--no-homology"] if ring else [])
    t = time.time()
    try:
        p = subprocess.run(args, capture_output=True, text=True, timeout=LIMIT)
    except subprocess.TimeoutExpired:
        return None, LIMIT
    if p.returncode: raise RuntimeError(p.stderr)
    return json.loads(p.stdout), time.time() - t


def hap_homology(g, D):
    out, t = gap('LoadPackage("hap");; SizeScreen([4096, 24]);; G := SmallGroup(%d, %d);; t := Runtime();;\n'
                 'R := ResolutionFiniteGroup(G, %d);; C := TensorWithIntegers(R);; H := List([1..%d], k -> Homology(C, k));;\n'
                 'Print("RESULT ", Runtime() - t, " ", H, "\\n");' % (g[0], g[1], D + 1, D))
    if out is None or "RESULT" not in out: return None, t, None
    m = re.search(r"RESULT (\d+) (.*)", out.replace("\\\n", "").replace("\n", " "))
    return json.loads(m.group(2).replace(" ", "")), t, int(m.group(1)) / 1000


def hap_ring(g, D):
    out, t = gap('LoadPackage("hap");; G := SmallGroup(%d, %d);; t := Runtime();;\n'
                 'P := Mod2CohomologyRingPresentation(G, %d);; Print("RESULT ", Runtime() - t, "\\n");' % (g[0], g[1], D))
    if out is None or "RESULT" not in out: return None, t, None
    return True, t, int(re.search(r"RESULT (\d+)", out).group(1)) / 1000


def main():
    quick = "--quick" in sys.argv
    outpath = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else os.path.join(HERE, "vs_hap.json")
    rows = []
    for g, D in (QUICK_H if quick else HOMOLOGY):
        f = export(g)
        pj, pt = perseus(f, D)
        hh, ht, hcpu = hap_homology(g, D)
        if pj is not None and hh is not None:
            ours = [h["torsion"] for h in pj["homology"][1:]]
            if ours != [sorted(x) for x in hh]: raise SystemExit("MISMATCH on SmallGroup%s: perseus %s, HAP %s" % (g, ours, hh))
        row = {"group": "SmallGroup(%d, %d)" % g, "task": "H_n(G; Z), n <= %d" % D, "perseus_s": round(pt, 3),
               "hap_s": None if hh is None else round(ht, 3), "hap_cpu_s": hcpu, "agree": pj is not None and hh is not None,
               "chains_top": pj["chains"][-1] if pj else None}
        rows.append(row); print(row, flush=True)
    for g, D in (QUICK_R if quick else RING):
        f = export(g)
        pj, pt = perseus(f, D, ring=2)
        hr, ht, hcpu = hap_ring(g, D)
        row = {"group": "SmallGroup(%d, %d)" % g, "task": "mod-2 ring, degrees <= %d" % D, "perseus_s": round(pt, 3),
               "hap_s": None if hr is None else round(ht, 3), "hap_cpu_s": hcpu,
               "perseus_generators_by_degree": pj["cohomology"]["generators_by_degree"][1:] if pj else None}
        rows.append(row); print(row, flush=True)
    json.dump({"limit_s": LIMIT, "rows": rows}, open(outpath, "w"), indent=1)


if __name__ == "__main__":
    main()
