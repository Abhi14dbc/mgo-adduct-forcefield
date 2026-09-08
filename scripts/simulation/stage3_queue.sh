#!/bin/bash
# STAGE 3 — extend the six apo-SH replicates from 100 ns to 300 ns.
#
# Why: at 100 ns the apo-SH replicates are not converged (RMSF half-to-half
# correlation r = 0.51-0.84, against 0.82-0.90 for apo-SS), so the apo-SH values
# in the manuscript are currently provisional. Extending to 300 ns also makes
# apo-SH directly comparable with apo-SS at the SAME window, so prediction P1
# can be tested at both a matched 100 ns and a matched 300 ns window instead of
# only the shorter one.
#
# Cost: 6 replicates x 200 ns = 1200 ns. At the measured 390 ns/day that is
# about 3.1 GPU-days.
#
# Ordering: this waits for stage 2 to finish, so the dimer condition (P2, P3 -
# the most disease-relevant result outstanding) is not delayed.
#
# Logic lives in this FILE, never `bash -c`, so polling for run_replicate.py
# cannot match the poller's own command line.
export HOME=/home/jyoti
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate mgo-md

ROOT=/home/jyoti/Projects/MGO
cd "$ROOT"
export MD_PLATFORM=CUDA
export PYTHONWARNINGS=ignore
mkdir -p logs
STEPS=150000000        # 300 ns

wait_for_stage2 () {
  while pgrep -f "bash stage2_queue.sh" > /dev/null || pgrep -f "run_replicate.py" > /dev/null; do
    sleep 120
  done
}

run () {   # run <name> <prmtop> <inpcrd> <seed>
  local name=$1 prm=$2 crd=$3 seed=$4
  local have=0
  [ -f "$name/md.log" ] && have=$(tail -1 "$name/md.log" | awk -F, '{printf "%.0f", $2*0.002/1000}')
  if [ "$have" -ge 299 ] 2>/dev/null; then
    echo "[$(date +%F_%T)] $name already at ${have} ns - skipping"; return
  fi
  echo "[$(date +%F_%T)] >>> $name (seed $seed) at ${have} ns -> 300 ns"
  python run_replicate.py "$prm" "$crd" "$name" "$seed" "$STEPS" \
      >> "logs/${name}_300.out" 2>&1
  echo "[$(date +%F_%T)] <<< $name exited $?"
}

echo "=== STAGE 3 queued $(date +%F_%T) - apo-SH 100 -> 300 ns, waiting for stage 2 ==="
wait_for_stage2
echo "=== STAGE 3 started $(date +%F_%T) ==="

SH=$ROOT/sod1_new
for i in 1 2 3; do
  run "nativeSH_rep$i" "$SH/native_SH.prmtop" "$SH/native_SH_final.inpcrd" "400$i"
done
for i in 1 2 3; do
  run "glycSH_rep$i"   "$SH/glyc_SH.prmtop"   "$SH/glyc_SH_final.inpcrd"   "500$i"
done

echo "=== STAGE 3 COMPLETE $(date +%F_%T) - apo-SH at 300 ns ==="
