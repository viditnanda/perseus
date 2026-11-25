#!/bin/sh
# Filtrations you build yourself: simplicial complexes and weighted graphs. Run:  sh examples/cli/04_your_own_filtration.sh
set -e
P=${PERSEUS:-perseus}
cd "$(dirname "$0")/../data"

echo "== lower-star filtration of a height function on a triangulated torus: one simplex per line, then its birth"
head -n 4 torus_height.txt
echo "   ..."
"$P" simplices torus_height.txt

echo
echo "== clique filtration of a weighted graph: lines 'u v birth' (edges) and 'v birth' (vertices); labels are words"
cat square_flag.txt
"$P" flag square_flag.txt --maxdim 2
