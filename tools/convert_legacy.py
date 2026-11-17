#!/usr/bin/env python3
"""Convert Perseus 4/5 input files to Perseus 6 formats, with exactly the old semantics.

    python3 convert_legacy.py <old format> <old file> <new file>

It writes the new file and prints the Perseus 6 command that computes the same filtration. Perseus 4/5
computed over F_2, so the command includes --coeff 2; drop it to get the rational barcode with its bad
primes.

  simtop, nmfsimtop (and _static, _capped)  ->  simplices
  cubtop                                    ->  cubical   (+ --periodic for negative extents)
  scubtop                                   ->  cubes
  cellcomp                                  ->  cells
  distmat, corrmat, brips, rips             ->  flag      (+ --maxdim: Perseus 4/5 built simplices up to a
                                                           dimension cap C; degrees below C are exact)
"""
import math, sys, re


def tokens(path):
    out = []
    for line in open(path):
        line = line.split("#", 1)[0]
        out.extend(t for t in re.split(r"[\s,]+", line) if t)
    return out


class Reader:
    def __init__(self, path): self.t = tokens(path); self.i = 0
    def more(self): return self.i < len(self.t)
    def num(self):
        if self.i >= len(self.t): raise EOFError
        self.i += 1; return float(self.t[self.i - 1])
    def int(self):
        x = self.num(); return int(x)


def fmt(x):
    return repr(int(x)) if float(x).is_integer() else repr(float(x))


def simtop(path, out, uniform, births, capped):
    r = Reader(path)
    topdim = r.int() if uniform else None
    ptdim = r.int()
    cap = r.int() if capped else None
    vid, simp = {}, {}
    try:
        while r.more():
            d = topdim if uniform else r.int()
            if d < 0: continue
            verts = []
            for _ in range(d + 1):
                pt = tuple(r.num() for _ in range(ptdim))
                verts.append(vid.setdefault(pt, len(vid)))
            b = r.int() if births else 0
            s = tuple(sorted(set(verts)))
            simp[s] = min(b, simp.get(s, b))          # duplicates: the oldest survives
    except EOFError:
        pass
    with open(out, "w") as f:
        for s, b in simp.items():
            if cap is not None and len(s) - 1 > cap:
                from itertools import combinations
                for t in combinations(s, cap + 1): f.write(" ".join(map(str, t)) + f" {b}\n")
            else:
                f.write(" ".join(map(str, s)) + f" {b}\n")
    return ["simplices", out, "--coeff", "2"]


def cubtop(path, out):
    r = Reader(path)
    d = r.int(); ext = []
    for _ in range(d):
        e = r.int()
        if e != 0: ext.append(e)
    shape = [abs(e) for e in ext]; tot = math.prod(shape)
    vals = []
    for _ in range(tot):
        try: vals.append(r.int())
        except EOFError: vals.append(-1)
    shift = 1 if 0 in vals else 0                     # Perseus: a birth of 0 shifts every birth by +1
    # old order: axis 0 fastest = C order of the reversed shape
    with open(out, "w") as f:
        f.write(" ".join(map(str, reversed(shape))) + "\n")
        f.write("\n".join("nan" if v == -1 else str(v + shift) for v in vals) + "\n")
    per = [str(len(ext) - 1 - k) for k, e in enumerate(ext) if e < 0]
    return ["cubical", out, "--coeff", "2"] + (["--periodic", ",".join(sorted(per))] if per else [])


def scubtop(path, out):
    r = Reader(path)
    d = r.int(); top = {}
    try:
        while r.more():
            a = tuple(r.int() for _ in range(d)); b = r.int()
            if b == -1: continue
            top.setdefault(a, b)                         # Perseus: the first occurrence of an anchor wins
    except EOFError:
        pass
    with open(out, "w") as f:
        for a, b in top.items(): f.write(" ".join(map(str, a)) + f" {b}\n")
    return ["cubes", out, "--coeff", "2"]


def cellcomp(path, out):
    r = Reader(path)
    ids, lines, cur = [], [], -1
    while r.more():
        x = r.int()
        if x < 0: cur += 1; ids.append([]); continue
        nb = r.int(); faces = []
        for _ in range(nb):
            c = r.int(); ix = r.int(); faces.append((ids[cur - 1][ix], c))
        ids[cur].append(len(lines))
        lines.append(f"{cur} {x}" + "".join(f" {g} {c}" for g, c in faces))
    open(out, "w").write("\n".join(lines) + "\n")
    return ["cells", out, "--coeff", "2"]


def write_flag(out, n, vb, edges):
    with open(out, "w") as f:
        for i in range(n):
            f.write(f"{i} {fmt(vb[i])}\n")              # every vertex, so isolated ones exist
        for (i, j), b in edges.items(): f.write(f"{i} {j} {fmt(b)}\n")


def distmat(path, out, corr):
    r = Reader(path)
    n = r.int(); g, s, N, d = r.num(), r.num(), r.int(), r.int()
    edges = {}
    for i in range(n):
        for j in range(n):
            x = r.num()
            if j >= i: continue
            if corr: x = 1 - x
            st = 0 if x < g else math.ceil((x - g) / s)
            if st <= N: edges[(i, j)] = st
    write_flag(out, n, [0] * n, edges)
    return ["flag", out, "--coeff", "2", "--maxdim", str(max(0, d - 1))]


def brips(path, out):
    r = Reader(path)
    pd = r.int(); g, s, N = r.num(), r.num(), r.int()
    pts = {}
    try:
        while r.more():
            p = tuple(r.num() for _ in range(pd)); rad = r.num()
            if rad >= 0: pts.setdefault(p, rad)
    except EOFError:
        pass
    P = sorted(pts); R = [pts[p] for p in P]; n = len(P)
    edges = {}
    for i in range(n):
        for j in range(i):
            st = max(0, math.ceil(0.5 / s * (math.dist(P[i], P[j]) / g - (R[i] + R[j]))))
            if st <= N: edges[(i, j)] = max(1, st)      # vertices are born at step 1
    write_flag(out, n, [1] * n, edges)
    return ["flag", out, "--coeff", "2", "--maxdim", str(max(0, pd - 1))]


def rips(path, out):
    r = Reader(path)
    pd = r.int(); pts, order = {}, []
    try:
        while r.more():
            p = tuple(r.num() for _ in range(pd)); rad = r.num()
            if rad < 0: continue
            b = r.num()
            if b == -1: continue
            if p not in pts: pts[p] = (rad, b); order.append(p)
    except EOFError:
        pass
    n = len(order); vb = [pts[p][1] for p in order]; edges = {}
    for i in range(n):
        for j in range(i):
            if math.dist(order[i], order[j]) < pts[order[i]][0] + pts[order[j]][0]:
                edges[(i, j)] = max(vb[i], vb[j])
    write_flag(out, n, vb, edges)
    return ["flag", out, "--coeff", "2", "--maxdim", str(max(0, pd - 1))]


def convert(fmt_, path, out):
    f = fmt_.lower()
    if f in ("simtop", "simtop_static", "nmfsimtop", "nmfsimtop_static", "nmfsimtop_capped", "nmfsimtop_static_capped"):
        return simtop(path, out, f.startswith("simtop"), "static" not in f, "capped" in f)
    if f == "cubtop": return cubtop(path, out)
    if f == "scubtop": return scubtop(path, out)
    if f == "cellcomp": return cellcomp(path, out)
    if f in ("distmat", "corrmat"): return distmat(path, out, f == "corrmat")
    if f == "brips": return brips(path, out)
    if f == "rips": return rips(path, out)
    raise SystemExit(f"unknown old format {fmt_}")


if __name__ == "__main__":
    if len(sys.argv) != 4: raise SystemExit(__doc__)
    args = convert(*sys.argv[1:])
    print("perseus " + " ".join(args))
