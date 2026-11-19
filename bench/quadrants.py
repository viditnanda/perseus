#!/usr/bin/env python3
"""Count value quadrants with Delta(s,t) != +-1 from a PERSEUS_DUMP file (k num den birth death).

Values: the sorted distinct filtration values v_0 < ... < v_N (file argument 2: one value per line).
Quadrant (s, t), -1 <= s < t <= N, of boundary degree k+1: rows = k-cells born after v_s, columns =
(k+1)-cells born by v_t. It contains the pair (b, d) iff ib > s and id <= t, where ib, id are value
indices. Delta != +-1 iff it contains a non-unit pair."""
import sys, bisect, collections

dump, valfile = sys.argv[1], sys.argv[2]
vals = [float(x) for x in open(valfile)]
N = len(vals) - 1
idx = lambda v: bisect.bisect_left(vals, v - 1e-7 * max(1.0, abs(v)))
pairs = collections.defaultdict(list)
for line in open(dump):
    k, num, den, b, d = line.split()
    d = float(d)
    if d == float("inf"):
        continue
    pairs[int(k)].append((idx(float(b)), idx(d)))
total = (N + 1) * (N + 2) // 2          # pairs -1 <= s < t <= N
out = {"values": N + 1, "quadrants_total": total}
for k, P in sorted(pairs.items()):
    # area of the union of {s < ib, t >= id}: for each s, T(s) = min id over pairs with ib > s
    P.sort()
    ibs = sorted({ib for ib, _ in P})
    # suffix minimum of id over pairs with ib > s
    bymin = {}
    for ib, idd in P:
        bymin[ib] = min(bymin.get(ib, 10**18), idd)
    area = 0
    suf = 10**18
    # iterate s from N-1 down to -1; pairs with ib > s
    keys = sorted(bymin, reverse=True)
    j = 0
    for s in range(N - 1, -2, -1):
        while j < len(keys) and keys[j] > s:
            suf = min(suf, bymin[keys[j]]); j += 1
        if suf <= N:
            area += N - max(suf, s + 1) + 1
    # distinct nonempty sets of non-unit pairs: sets change only when s crosses an ib or t an id
    ids = sorted({idd for _, idd in P})
    seen = set()
    for ib0 in ibs:
        s = ib0 - 1
        for t in ids:
            S = frozenset(i for i, (ib, idd) in enumerate(P) if ib > s and idd <= t)
            if S:
                seen.add(S)
    out[f"deg{k}"] = dict(nonunit=len(P), quadrants_nonunit=area, fraction=area / total,
                          distinct_sets=len(seen), row_thresholds=len(ibs), col_thresholds=len(ids))
import json
print(json.dumps(out))
