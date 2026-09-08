#!/bin/bash
# STAGE 2 - the three supporting studies, after the six-run 300 ns set.
#
# Order is deliberate: the two cheap, high-value studies run first so the paper
# can be submitted without waiting on the expensive one.
#
#   1. apo-SH        1.6 d  removes the "wrong species" objection - apo-SH is
#                           both the most aggregation-prone form (Toichi 2013)
#                           and the most MGO-reactive (Polykretis 2020)
#   2. adduct series 1.6 d  MG-H1 vs argpyrimidine vs CEA at ONE site, charge
#                           held constant - demonstrates why parameterising four
#                           adducts mattered; no prior study exists
#   3. dimer         3.4 d  biggest scientific claim, highest cost - last, so it
#                           gates nothing
#
# New conditions run to 100 ns. The four-way adduct comparison then uses a
# matched 10-100 ns window across native/MG-H1/ARP/CEA; the 300 ns data serves
# the primary two-way native-vs-glycated comparison.
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
STEPS=50000000        # 100 ns

wait_for_stage1 () {
  # stage 1 is the master queue; wait for it AND any live replicate
  while pgrep -f "bash master_queue.sh" > /dev/null || pgrep -f "run_replicate.py" > /dev/null; do
    sleep 120
  done
}

run () {   # run <name> <prmtop> <inpcrd> <seed>
  local name=$1 prm=$2 crd=$3 seed=$4
  local have=0
  [ -f "$name/md.log" ] && have=$(tail -1 "$name/md.log" | awk -F, '{printf "%.0f", $2*0.002/1000}')
  if [ "$have" -ge 99 ] 2>/dev/null; then
    echo "[$(date +%F_%T)] $name already at ${have} ns - skipping"; return
  fi
  echo "[$(date +%F_%T)] >>> $name (seed $seed)"
  python run_replicate.py "$prm" "$crd" "$name" "$seed" "$STEPS" \
      >> "logs/${name}.out" 2>&1
  echo "[$(date +%F_%T)] <<< $name exited $?"
}

echo "=== STAGE 2 queued $(date +%F_%T) - waiting for stage 1 ==="
wait_for_stage1
echo "=== STAGE 2 started $(date +%F_%T) ==="

SH=$ROOT/sod1_new
AD=$ROOT/sod1_adducts

# ---- 1. apo-SH (disulfide reduced) ----
for i in 1 2 3; do
  run "nativeSH_rep$i" "$SH/native_SH.prmtop" "$SH/native_SH_final.inpcrd" "400$i"
done
for i in 1 2 3; do
  run "glycSH_rep$i"   "$SH/glyc_SH.prmtop"   "$SH/glyc_SH_final.inpcrd"   "500$i"
done

# ---- 2. adduct series at Arg143 ----
for i in 1 2 3; do
  run "arp143_rep$i" "$AD/glyc_SOD1apo_ARP143.prmtop" "$AD/glyc_SOD1apo_ARP143_L.inpcrd" "600$i"
done
for i in 1 2 3; do
  run "cea143_rep$i" "$AD/glyc_SOD1apo_CEA143.prmtop" "$AD/glyc_SOD1apo_CEA143_L.inpcrd" "700$i"
done

# ---- 2b. CEL at the ALS-implicated lysines ----
# Monteiro Neto et al. (2023) name Lys122/Lys128 in sporadic ALS models. CEL
# accumulates in long-lived proteins whereas MG-H1 (t-half ~12 d) does not, so
# for a chronic disease this is arguably the more relevant adduct class. CEL
# also changes charge by -2 per site (Lys+ -> carboxylate-) versus -1 for the
# arginine adducts.
CEL=$ROOT/sod1_cel
for i in 1 2 3; do
  run "cel122_rep$i" "$CEL/cel122.prmtop" "$CEL/cel122_L.inpcrd" "710$i"
done
for i in 1 2 3; do
  run "cel128_rep$i" "$CEL/cel128.prmtop" "$CEL/cel128_L.inpcrd" "720$i"
done

# ---- 3. dimer ----
for i in 1 2 3; do
  run "nativeDim_rep$i" "$SH/native_dimer.prmtop" "$SH/native_dimer_final.inpcrd" "800$i"
done
for i in 1 2 3; do
  run "glycDim_rep$i"   "$SH/glyc_dimer.prmtop"   "$SH/glyc_dimer_final.inpcrd"   "900$i"
done

echo "=== STAGE 2 COMPLETE $(date +%F_%T) ==="
