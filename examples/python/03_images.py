"""Cubical persistence of images and volumes, straight from numpy.

Each array entry is the value of one pixel (voxel); the filtration adds pixels in increasing order of value.
NaN marks an absent pixel. Run:  python3 examples/python/03_images.py
"""
import os

import numpy as np
import perseus

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


def long_bars(barcode, at_least):
    """The intervals of each degree whose length is at least `at_least` (essential ones included)."""
    return {k: np.round(B[B[:, 1] - B[:, 0] >= at_least], 3).tolist() for k, B in enumerate(barcode) if len(B)}


# A 64 x 64 image: two dark discs and a dark ring on a bright background (examples/data/make_data.py).
img = np.load(os.path.join(DATA, "blobs.npy"))
r = perseus.cubical(img)
print("image", img.shape, "| certificate:", r["certificate"]["arithmetic"])
print("  all intervals:", [len(B) for B in r["barcode"]], "(H0, H1)")
print("  intervals of length >= 0.1:", long_bars(r["barcode"], 0.1))
print("  -> three components (two discs and the ring; two of them merge into the third at the background level)")
print("     and one loop (the hole of the ring)")

# The superlevel filtration (bright features first): negate the image.
r = perseus.cubical(-img)
print("  superlevel sets, intervals of length >= 0.1:", long_bars(r["barcode"], 0.1))

# Periodic boundary conditions. Periodic in one axis, the grid is a cylinder, and a stripe that wraps
# around it is a loop.
stripe = np.ones((32, 32))
stripe[10:14, :] = 0.0
r = perseus.cubical(stripe, periodic=[False, True])
print("\nstripe, periodic in axis 1:", long_bars(r["barcode"], 0.5))

# Absent pixels: NaN (or inf). The essential classes of what remains are reported as infinite intervals.
frame = np.zeros((10, 10))
frame[3:7, 3:7] = np.nan
print("a frame with a hole of absent pixels:", long_bars(perseus.cubical(frame)["barcode"], 0))

# A 3D volume. In 2D the planar engine runs; in 3D and 4D the discrete-gradient (Morse) engine.
vol = np.load(os.path.join(DATA, "volume.npy"))
r = perseus.cubical(vol)
print("\nvolume", vol.shape, ": intervals per degree", [len(B) for B in r["barcode"]])
print("  H0 intervals of length >= 0.2:", len(long_bars(r["barcode"], 0.2)[0]), "(the deepest blobs; the others merge early)")
print("  bad primes:", perseus.bad_primes(r) or "none")

# Engines can be chosen by hand, e.g. to compare: auto | planar | morse | classic | explicit.
same = all(np.array_equal(np.sort(a, axis=0), np.sort(b, axis=0))
           for a, b in zip(perseus.cubical(vol, engine="classic", coeff=2)["barcode"],
                           perseus.cubical(vol, engine="morse", coeff=2)["barcode"]))
print("  classic and Morse engines agree over F_2:", same)
