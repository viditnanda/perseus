#!/usr/bin/env python3
"""Regenerate the generated example inputs in this directory (numpy only, fixed seeds).

    python3 examples/data/make_data.py

The hand-written inputs (image_cubical.txt, moore12_cells.txt, rp2_simplices.txt, square_flag.txt) and the
finite-group files (written from GAP by tools/export_*.g) are not touched.
"""
import os
from itertools import combinations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def path(name):
    return os.path.join(HERE, name)


def circle():
    """100 points near the unit circle: one long H1 bar."""
    rng = np.random.default_rng(1)
    t = rng.uniform(0, 2 * np.pi, 100)
    X = np.c_[np.cos(t), np.sin(t)] + rng.normal(scale=0.08, size=(100, 2))
    with open(path("circle.txt"), "w") as fh:
        fh.write("# 100 points near the unit circle (examples/data/make_data.py)\n")
        np.savetxt(fh, X, fmt="%.6f")


def rp2_points():
    """RP^2 as a point cloud: the Veronese image x -> (x_i x_j) of 250 points of the 2-sphere, in R^6.
    Antipodal points have the same image, so this samples the real projective plane."""
    rng = np.random.default_rng(2)
    S = rng.normal(size=(250, 3))
    S /= np.linalg.norm(S, axis=1)[:, None]
    V = np.c_[S[:, 0] ** 2, S[:, 1] ** 2, S[:, 2] ** 2,
              np.sqrt(2) * S[:, 0] * S[:, 1], np.sqrt(2) * S[:, 0] * S[:, 2], np.sqrt(2) * S[:, 1] * S[:, 2]]
    with open(path("rp2_points.txt"), "w") as fh:
        fh.write("# 250 points on the Veronese RP^2 in R^6 (examples/data/make_data.py)\n")
        np.savetxt(fh, V, fmt="%.6f")


def images():
    """A 64 x 64 grey-scale image (two dark discs and a dark ring on a bright background; three components
    and one hole in its sublevel sets) and a 32^3 volume (a dozen dark blobs and a little noise)."""
    y, x = np.mgrid[0:64, 0:64].astype(float)
    bump = lambda d2, s: np.exp(-d2 / (2 * s * s))
    img = 1 - np.maximum.reduce([
        bump((x - 14) ** 2 + (y - 14) ** 2, 4),          # a disc
        bump((x - 14) ** 2 + (y - 50) ** 2, 3),          # a smaller disc
        bump((np.hypot(x - 40, y - 40) - 12) ** 2, 2),   # a ring around (40, 40)
    ]) + 1e-3 * (x + 0.37 * y) / 64                      # a slight tilt, so that no two pixels tie
    np.save(path("blobs.npy"), img)
    rng = np.random.default_rng(3)
    centres = rng.uniform(4, 28, size=(12, 3))
    z, y, x = np.mgrid[0:32, 0:32, 0:32].astype(float)
    vol = np.zeros((32, 32, 32))
    for c in centres:
        vol -= np.exp(-((x - c[0]) ** 2 + (y - c[1]) ** 2 + (z - c[2]) ** 2) / 18.0)
    vol += 0.02 * rng.standard_normal(vol.shape)
    np.save(path("volume.npy"), vol)


def torus_lower_star():
    """A triangulated torus (an 8 x 8 grid with opposite sides identified) filtered by height:
    each vertex is born at its height in the standard embedding, each simplex at its highest vertex."""
    n = 8
    def v(i, j):
        return (i % n) * n + (j % n)
    height = {}
    for i in range(n):
        for j in range(n):
            a, b = 2 * np.pi * i / n, 2 * np.pi * j / n
            height[v(i, j)] = round(float((2 + np.cos(b)) * np.cos(a)), 6)   # x-coordinate: a Morse function
    tris = []
    for i in range(n):
        for j in range(n):
            tris.append((v(i, j), v(i + 1, j), v(i + 1, j + 1)))
            tris.append((v(i, j), v(i, j + 1), v(i + 1, j + 1)))
    edges = sorted({tuple(sorted(e)) for t in tris for e in combinations(t, 2)})
    with open(path("torus_height.txt"), "w") as fh:
        fh.write("# torus: 8 x 8 grid triangulation, lower-star filtration of a height function\n")
        fh.write("# each line: the vertices of a simplex, then its birth (the largest height of its vertices)\n")
        for k, h in sorted(height.items()):
            fh.write(f"{k} {h}\n")
        for s in edges + tris:
            fh.write(" ".join(map(str, s)) + f" {max(height[u] for u in s)}\n")


def rp2_two_stage():
    """RP^1 at time 0, then RP^2 at time 1, as a based chain complex over Z (one cell per dimension)."""
    with open(path("rp2_two_stage.cells"), "w") as fh:
        fh.write("# K_0 = RP^1 (a vertex and a loop), K_1 = RP^2 (a 2-cell attached by degree 2)\n")
        fh.write("# each line: dim birth [face coefficient]...; faces are 0-based line numbers of earlier cells\n")
        fh.write("0 0\n1 0\n2 1 1 2\n")


def linial_meshulam():
    """A Linial-Meshulam process on 40 vertices: all edges at time 1, then 700 random triangles one at a time.
    Near the end H_1 dies, and on its way out it carries torsion (Z/975 = Z/(3 * 5^2 * 13) at one point)."""
    rng = np.random.default_rng(40)
    n, m = 40, 700
    tris = list(combinations(range(n), 3))
    order = rng.permutation(len(tris))[:m]
    with open(path("lm40.txt"), "w") as fh:
        fh.write("# Linial-Meshulam process: all edges of K_40 at time 1, then 700 random triangles at times 2, 3, ...\n")
        for a, b in combinations(range(n), 2):
            fh.write(f"{a} {b} 1\n")
        for t, i in enumerate(order):
            fh.write(" ".join(map(str, tris[i])) + f" {t + 2}\n")


def simplicial_set():
    """RP^2 as a simplicial set with one nondegenerate simplex in each dimension 0, 1, 2."""
    with open(path("rp2.sset"), "w") as fh:
        fh.write("# RP^2 as a simplicial set: vertex v, edge e (d0 e = d1 e = v), triangle f with faces e, s_0 v, e\n")
        fh.write("# each line: dim birth d_0 ... d_dim; a face is an earlier line, y@k for the degenerate s_k y\n")
        fh.write("0 0\n1 0 0 0\n2 1 1 0@0 1\n")


def magnitude_graphs():
    """The 7-cycle, and the Kaneta-Yoshinaga graph of the 6-vertex RP^2 (the Hasse diagram of its face poset
    with a bottom and a top added), whose magnitude homology has 2-torsion."""
    with open(path("cycle7.txt"), "w") as fh:
        fh.write("# the cycle graph C_7: one edge u v per line, vertices 0..6\n")
        for i in range(7):
            fh.write(f"{i} {(i + 1) % 7}\n")
    facets = [(0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5), (0, 1, 5), (1, 2, 4), (2, 3, 5), (1, 3, 4), (2, 4, 5), (1, 3, 5)]
    faces = sorted({tuple(sorted(s)) for F in facets for k in (1, 2, 3) for s in combinations(F, k)}, key=lambda s: (len(s), s))
    idx = {s: i for i, s in enumerate(faces)}
    bot, top = len(faces), len(faces) + 1
    E = []
    for s in faces:
        if len(s) == 1:
            E.append((bot, idx[s]))
        if len(s) == 3:
            E.append((idx[s], top))
        else:
            for u in faces:
                if len(u) == len(s) + 1 and set(s) <= set(u):
                    E.append((idx[s], idx[u]))
    with open(path("ky_rp2.txt"), "w") as fh:
        fh.write(f"# Kaneta-Yoshinaga graph of the 6-vertex RP^2: {len(faces)} faces, bottom {bot}, top {top}\n")
        for a, b in E:
            fh.write(f"{a} {b}\n")


def trefoil():
    """The infinite cyclic cover of the trefoil complement, from the presentation <x, y | xyx = yxy>
    by Fox calculus: a chain complex over Z[t, t^-1] whose H_1 is Z[t, t^-1]/(t^2 - t + 1)."""
    with open(path("trefoil.laurent"), "w") as fh:
        fh.write("# trefoil complement, infinite cyclic cover: one vertex, generators x, y, relator xyx(yxy)^-1\n")
        fh.write("# each line: dim birth face:polynomial ...; faces are 0-based line numbers of earlier cells\n")
        fh.write("0 0\n1 0 0:t-1\n1 0 0:t-1\n2 0 1:t^2-t+1 2:-t^2+t-1\n")


if __name__ == "__main__":
    circle(); rp2_points(); images(); torus_lower_star(); rp2_two_stage(); linial_meshulam()
    simplicial_set(); magnitude_graphs(); trefoil()
