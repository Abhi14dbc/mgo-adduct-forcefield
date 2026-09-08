#!/bin/bash
# STAGE 4 — take the apo-SH implicit arm from n = 3 to n = 10.
#
# Why: the solvent-model comparison is currently asymmetric. The explicit-OPC
# arm on KISTI has n = 10 (permutation floor 1.1e-05); this implicit-GB arm has
# n = 3 (floor 0.10), so no endpoint here can reach significance at any alpha
# no matter how large the effect. That is a property of the design, not of the
# data, and it is the one thing a reviewer can attack the comparison on.
#
# The two arms are otherwise identical and verified so: both apo (no Zn/Cu),
# both disulfide-reduced (CYS at 5/56/110/145, no CYX), both 153 protein
# residues, same adduct at index 142, analysed by the same script at the same
# 200 ps frame resolution over the same final-150 ns window. Solvent model is
# the ONLY variable. Bringing this arm to n = 10 makes it symmetric.
#
# Protocol is UNCHANGED from stage 2/3 — same integrator, thermostat, cutoff,
# timestep, starting coordinates (*_final.inpcrd) and seed formula. Replicates
# 4-10 differ from 1-3 only in velocity seed, which is what makes all ten
# poolable as one sample. Do not "improve" anything here.
#
# Seeds: 4000+i native, 5000+i glycated. Reps 1-3 used 4001-4003 / 5001-5003
# under the same formula, so 4004-4010 / 5004-5010 continue it without collision.
#
# Cost: 14 replicates x 300 ns = 4200 ns. At the measured ~380 ns/day that is
# about 11 GPU-days on the local card. Zero KISTI SRU.
#
# Idempotent: a replicate already at >= 299 ns is skipped, and run_replicate.py
# auto-resumes from state.chk, so this survives interruption and can be re-run.
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
STEPS=150000000        # 300 ns at 2 fs

wait_for_running () {
  while pgrep -f "run_replicate.py" > /dev/null; do
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

echo "=== STAGE 4 queued $(date +%F_%T) - apo-SH n=3 -> n=10, waiting for any live run ==="
wait_for_running
echo "=== STAGE 4 started $(date +%F_%T) ==="

SH=$ROOT/sod1_new
for i in 4 5 6 7 8 9 10; do
  run "nativeSH_rep$i" "$SH/native_SH.prmtop" "$SH/native_SH_final.inpcrd" "$((4000+i))"
done
for i in 4 5 6 7 8 9 10; do
  run "glycSH_rep$i"   "$SH/glyc_SH.prmtop"   "$SH/glyc_SH_final.inpcrd"   "$((5000+i))"
done

echo "=== STAGE 4 COMPLETE $(date +%F_%T) - apo-SH at n=10, 300 ns ==="
echo "Next: edit GROUPS in analyse_stage3.py to range(1, 11) and re-run:"
echo "  /home/jyoti/miniforge3/envs/mgo-md/bin/python analyse_stage3.py"
