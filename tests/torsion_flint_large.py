#!/usr/bin/env python3
"""Independent check of `perseus --torsion` on large Linial-Meshulam instances at chosen t (needs python-flint):
exact sparse elimination of unit entries over Z (unimodular, so the cokernel is unchanged), then a FLINT
Smith form of the small residual."""
import sys, json, time
import flint

path, ts, jpath = sys.argv[1], [int(x) for x in sys.argv[2].split(",")], sys.argv[3]
lines = [l.split() for l in open(path)]
edges = sorted(tuple(map(int, l[:2])) for l in lines if len(l) == 3)
ei = {e: i for i, e in enumerate(edges)}
tris = sorted(((float(l[3]), tuple(map(int, l[:3]))) for l in lines if len(l) == 4))
steps = json.load(open(jpath))["torsion"]["degrees"][0]["steps"]
def pred(t):
    cur = []
    for s in steps:
        if s["value"] <= t: cur = [int(x) for x in s["invariant_factors"]]
    return cur

def torsion(t):
    cols = {}
    rows = {}
    for j, (v, (a, b, c)) in enumerate(x for x in tris if x[0] <= t):
        col = {ei[(b, c)]: 1, ei[(a, c)]: -1, ei[(a, b)]: 1}
        cols[j] = col
        for r in col: rows.setdefault(r, set()).add(j)
    while True:
        best = None
        for j, col in cols.items():
            for r, v in col.items():
                if v in (1, -1):
                    cost = (len(col) - 1) * (len(rows[r]) - 1)
                    if best is None or cost < best[0]: best = (cost, j, r)
            if best and best[0] == 0: break
        if best is None: break
        _, j, r = best
        pc = cols.pop(j); a = pc[r]
        for k in list(rows[r]):
            if k == j: continue
            ck = cols[k]; f = ck[r] * a          # a = +-1, so ck[r] / a = ck[r] * a
            for rr, vv in pc.items():
                nv = ck.get(rr, 0) - f * vv
                if nv:
                    if rr not in ck: rows.setdefault(rr, set()).add(k)
                    ck[rr] = nv
                else:
                    if rr in ck: del ck[rr]; rows[rr].discard(k)
        for rr in pc: rows[rr].discard(j)
        del rows[r]
        for k in cols:
            cols[k].pop(r, None)
    live_rows = sorted(r for r in rows if rows[r])
    ri = {r: i for i, r in enumerate(live_rows)}
    cl = [c for c in cols.values() if c]
    if not cl or not live_rows: return []
    M = flint.fmpz_mat(len(live_rows), len(cl))
    for j, c in enumerate(cl):
        for r, v in c.items(): M[ri[r], j] = v
    S = M.snf()
    return sorted(abs(int(S[i, i])) for i in range(min(S.nrows(), S.ncols())) if abs(int(S[i, i])) > 1), (len(live_rows), len(cl))

bad = 0
for t in ts:
    t0 = time.time()
    got = torsion(t)
    tors, shape = (got if got else ([], (0, 0)))
    ok = tors == pred(t); bad += not ok
    print(f"t={t}: residual {shape[0]}x{shape[1]}, smith {tors} vs perseus {pred(t)} {'OK' if ok else 'MISMATCH'} ({time.time()-t0:.1f}s)", flush=True)
sys.exit(1 if bad else 0)
