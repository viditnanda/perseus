#!/bin/sh
# Chain complexes over Z[t, t^-1]: the trefoil's infinite cyclic cover and its Alexander polynomial.
# Run:  sh examples/cli/08_laurent.sh
set -e
P=${PERSEUS:-perseus}
cd "$(dirname "$0")/../data"

echo "== each line: dim birth face:polynomial ..."
cat trefoil.laurent
"$P" laurent trefoil.laurent
