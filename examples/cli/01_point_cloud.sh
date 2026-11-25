#!/bin/sh
# Vietoris-Rips persistence of a point cloud: 100 noisy points near a circle.
# Run from anywhere:  sh examples/cli/01_point_cloud.sh      (PERSEUS=/path/to/perseus if it is not on PATH)
set -e
P=${PERSEUS:-perseus}
cd "$(dirname "$0")/../data"
OUT=${TMPDIR:-/tmp}/perseus_examples; mkdir -p "$OUT"

echo "== the barcode on the screen (Ripser's format), with the certificate summary on stderr"
"$P" points circle.txt --maxdim 1 | head -n 8
echo "   ..."

echo "== the same into files: circle_0.txt, circle_1.txt (birth death, -1 = never dies), circle_betti.txt, circle_cert.json"
"$P" points circle.txt "$OUT/circle" --maxdim 1 --print none
echo "   H1 intervals longer than 0.5:"
awk '$2 == -1 || $2 - $1 > 0.5' "$OUT/circle_1.txt"

echo "== only over F_2 (no certificate), ignoring edges longer than 1.0"
"$P" points circle.txt --maxdim 1 --coeff 2 --threshold 1.0 --print perseus | awk '$1 == 1 && ($3 == -1 || $3 - $2 > 0.5)'
