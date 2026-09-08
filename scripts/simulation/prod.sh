#!/bin/bash
#SBATCH --job-name=mgo_prod
#SBATCH --partition=cas_v100nv_8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:4
#SBATCH --time=24:00:00
#SBATCH --comment="field=bio;appl=openmm"
#SBATCH -o prod_%j.out
#SBATCH -e prod_%j.err
# Usage: sbatch prod.sh native_mono:1 native_mono:2 native_mono:3 native_mono:4
# One replicate per GPU; packing measured at 99.4% of ideal on this queue.
# Fewer than 4 replicates: override so you are not billed for idle GPUs, e.g.
#   sbatch --gres=gpu:2 --cpus-per-task=8 prod.sh native_mono:5 glyc_mono:5
source /scratch/a2148a01/miniforge3/etc/profile.d/conda.sh
conda activate mgo
cd /scratch/a2148a01
echo "host=$(hostname)  start=$(date)  args=$*"
nvidia-smi --query-gpu=index,name --format=csv,noheader
echo
i=0
for SPEC in "$@"; do
    SYS="${SPEC%%:*}"
    REP="${SPEC##*:}"
    CUDA_VISIBLE_DEVICES=$i python prod.py --system "$SYS" --rep "$REP" \
        --ns 300 > "prod_${SYS}_rep${REP}.txt" 2>&1 &
    i=$((i+1))
done
wait
echo
grep -h -E "START|DONE" prod_*.txt
echo "finished $(date)"
