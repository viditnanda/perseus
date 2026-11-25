#!/bin/sh
# Coefficients matter: the barcode over Q, every bad prime p (where F_p gives a different barcode), and
# what changes. Run:  sh examples/cli/02_bad_primes.sh
set -e
P=${PERSEUS:-perseus}
cd "$(dirname "$0")/../data"
OUT=${TMPDIR:-/tmp}/perseus_examples; mkdir -p "$OUT"

echo "== the Moore space M(Z/12, 1): the loop dies over Q, but not over F_2 or F_3 (--fp prints those barcodes too)"
cat moore12_cells.txt
"$P" cells moore12_cells.txt --fp 2,3,5

echo
echo "== RP^2 (6 vertices): a point over Q, but H1 and H2 over F_2"
"$P" simplices rp2_simplices.txt --fp 2

echo
echo "== RP^2 as a cloud of 250 points in R^6; the certificate (JSON) says which intervals change over F_2"
"$P" points rp2_points.txt "$OUT/rp2" --maxdim 2 --print none
python3 -c "
import json, sys
c = json.load(open(sys.argv[1]))
print('proof:', c['proof'], '| bad primes:', [b['p'] for b in c['bad_primes']], '| complete:', c['complete'])
for b in c['bad_primes']:
    print('  only over F_%d:' % b['p'], b['only_over_Fp'])
    print('  only over Q:   ', b['only_over_Q'])
" "$OUT/rp2_cert.json" 2>/dev/null || grep -o '"bad_primes": \[[^]]*' "$OUT/rp2_cert.json"
