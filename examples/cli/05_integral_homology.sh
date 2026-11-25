#!/bin/sh
# Integral homology, images of maps, and the torsion along a filtration. Run:  sh examples/cli/05_integral_homology.sh
set -e
P=${PERSEUS:-perseus}
cd "$(dirname "$0")/../data"

echo "== RP^1 at time 0, RP^2 at time 1: H_*(K; Z) and the image of H_*(K_0; Z) -> H_*(K_1; Z)"
cat rp2_two_stage.cells
"$P" cells rp2_two_stage.cells --integral --image 0,1 --print none

echo
echo "== a Linial-Meshulam process on 40 vertices: the torsion of H_1(K_t; Z) at every t where it changes"
"$P" simplices lm40.txt --torsion --print none
