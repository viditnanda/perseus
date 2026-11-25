#!/bin/sh
# Simplicial sets: RP^2 with one nondegenerate simplex in each dimension. Run:  sh examples/cli/06_simplicial_sets.sh
# (Nerves of categories and bar constructions are in the Python module: examples/python/06_simplicial_sets.py.)
set -e
P=${PERSEUS:-perseus}
cd "$(dirname "$0")/../data"

echo "== each line: dim birth d_0 ... d_dim; a face is an earlier line, or y@k for the degenerate simplex s_k y"
cat rp2.sset
"$P" sset rp2.sset --integral --image 0,1
