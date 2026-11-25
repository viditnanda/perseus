#!/bin/sh
# Cubical persistence of images and volumes. Run:  sh examples/cli/03_images.sh
set -e
P=${PERSEUS:-perseus}
cd "$(dirname "$0")/../data"
OUT=${TMPDIR:-/tmp}/perseus_examples; mkdir -p "$OUT"

echo "== a 3 x 4 image as text: sizes on the first line, then the values (nan = absent pixel)"
cat image_cubical.txt
"$P" cubical image_cubical.txt

echo
echo "== a 64 x 64 image from numpy (blobs.npy: two dark discs and a dark ring); intervals longer than 0.1"
"$P" npy blobs.npy "$OUT/blobs" --print none
for k in 0 1; do echo "   H$k:"; awk '$2 == -1 || $2 - $1 > 0.1' "$OUT/blobs_$k.txt"; done

echo "== the same image, periodic in axis 1 (a cylinder): the background now closes up into a loop"
"$P" npy blobs.npy "$OUT/blobs_cyl" --periodic 1 --print none
echo "   H1:"; awk '$2 == -1 || $2 - $1 > 0.1' "$OUT/blobs_cyl_1.txt"

echo
echo "== a 32^3 volume: the discrete-gradient engine runs by default in 3D and 4D"
"$P" npy volume.npy "$OUT/volume" --print none
wc -l "$OUT"/volume_[0-9].txt | sed '$d'
