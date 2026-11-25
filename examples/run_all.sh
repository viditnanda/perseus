#!/bin/sh
# Run every example, command line and Python, and stop at the first failure. Used by `make examples` and CI.
#   PERSEUS=/path/to/perseus sh examples/run_all.sh          (the Python module must be importable)
HERE="$(cd "$(dirname "$0")" && pwd)"
LOG=${TMPDIR:-/tmp}/perseus_examples.log
for f in "$HERE"/cli/*.sh; do
    printf '%s ... ' "cli/$(basename "$f")"
    sh "$f" > "$LOG" 2>&1 || { echo FAILED; cat "$LOG"; exit 1; }
    echo ok
done
for f in "$HERE"/python/*.py; do
    printf '%s ... ' "python/$(basename "$f")"
    python3 "$f" > "$LOG" 2>&1 || { echo FAILED; cat "$LOG"; exit 1; }
    echo ok
done
echo "all examples ran"
