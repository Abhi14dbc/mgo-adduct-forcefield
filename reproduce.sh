#!/usr/bin/env bash
# Reproduce the reported values from the deposited trajectories.
#
# Usage:
#   ./reproduce.sh /path/to/trajectory/deposit
#
# The trajectory deposit is a separate record; see the paper's Data and code
# availability section for its DOI. Both analyses read their input directory from
# MGO_DIR, so nothing here copies or rewrites data.
#
# The deposited trajectories are at 200 ps per frame, which is the resolution at
# which every reported value was computed. The implicit-solvent script strides its
# input by 20 by default because the working trajectories are at 10 ps; on the
# deposit that stride is 1 and the input spacing is 200 ps. Their product -- the
# analysis frame spacing -- is 200 ps either way, which matters because switch
# rates scale with it.

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "usage: $0 /path/to/trajectory/deposit" >&2
    exit 2
fi

TRAJ="$(cd "$1" && pwd)"
PY="${PYTHON:-python3}"
HERE="$(cd "$(dirname "$0")" && pwd)"

for d in "$TRAJ/explicit_mono" "$TRAJ/nativeSH" "$TRAJ/glycSH"; do
    [ -d "$d" ] || { echo "missing $d -- is that the trajectory deposit?" >&2; exit 1; }
done

echo "=============================================================================="
echo "EXPLICIT SOLVENT (primary arm), n = 10 per group, final 150 ns"
echo "Expected: CORE composite native 65.94, glycated 29.80, delta -36.14, p = 0.0009"
echo "=============================================================================="
MGO_DIR="$TRAJ/explicit_mono" "$PY" "$HERE/scripts/analysis/analyse_phase1.py"

echo
echo "=============================================================================="
echo "IMPLICIT SOLVENT, n = 10 per group, final 150 ns"
echo "Expected: CORE composite native 47.89, glycated 58.15, delta +10.25, p = 0.4638"
echo "          site H-bonds delta -0.654, 95% CI -1.195 to -0.120, p = 0.0381"
echo "=============================================================================="
MGO_DIR="$TRAJ" MGO_RAW_PS=200 MGO_STRIDE=1 "$PY" "$HERE/scripts/analysis/analyse_stage3.py"

echo
echo "Both arms reproduced. Values above should match the paper exactly; the"
echo "bootstrap intervals are seeded (20260904) and are reproducible to the digit."
