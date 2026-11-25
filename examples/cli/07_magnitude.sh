#!/bin/sh
# Magnitude homology of graphs over Z, with torsion. Run:  sh examples/cli/07_magnitude.sh
set -e
P=${PERSEUS:-perseus}
cd "$(dirname "$0")/../data"

echo "== the 7-cycle (one edge 'u v' per line), lengths l <= 5"
"$P" magnitude cycle7.txt --lmax 5

echo
echo "== the Kaneta-Yoshinaga graph of RP^2 between its bottom (31) and top (32): torsion from H_1(RP^2) = Z/2"
"$P" magnitude ky_rp2.txt --lmax 6 --pairs 31,32
