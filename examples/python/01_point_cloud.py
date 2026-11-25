"""Vietoris-Rips persistence of a point cloud.

100 noisy points near a circle: the barcode in degree 1 has one long interval (the circle) and many
short ones (noise). Run from anywhere:  python3 examples/python/01_point_cloud.py
"""
import os

import numpy as np
import perseus

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")

# The same points as examples/data/circle.txt; any (n, d) array works.
X = np.loadtxt(os.path.join(DATA, "circle.txt"))

r = perseus.rips(points=X, maxdim=1)      # over Q, with the certificate of bad primes (here: none)
H0, H1 = r["barcode"][0], r["barcode"][1]  # (m, 2) arrays of [birth, death); death = inf if essential
print(f"{len(X)} points: {len(H0)} intervals in H0, {len(H1)} in H1")

lengths = H1[:, 1] - H1[:, 0]
order = np.argsort(-lengths)
print("the three longest H1 intervals:")
for b, d in H1[order[:3]]:
    print(f"  [{b:.3f}, {d:.3f})   length {d - b:.3f}")

# The certificate: is the barcode proved, and which primes p would give a different barcode over F_p?
c = r["certificate"]
print("proved:", c["proof"], "| bad primes:", perseus.bad_primes(r) or "none")

# A distance matrix works just as well (square, or condensed as from scipy.spatial.distance.pdist).
D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=-1)
r2 = perseus.rips(distances=D, maxdim=1, threshold=1.0)   # ignore edges longer than 1.0
H1t = r2["barcode"][1]
print("from the distance matrix, threshold 1.0: the circle is now", H1t[np.argmax(H1t[:, 1] - H1t[:, 0])].tolist(),
      "(it never dies below the threshold)")

# Over a single prime field, without the certificate: coeff=p.
r3 = perseus.rips(points=X, maxdim=1, coeff=2)
print("over F_2:", len(r3["barcode"][1]), "intervals in H1 (the same as over Q: there are no bad primes)")
