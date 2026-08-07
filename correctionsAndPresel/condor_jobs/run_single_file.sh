#!/bin/bash
set -e

source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc15-opt/setup.sh
cd /eos/home-a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/correctionsAndPresel

INPUT_FILE=$1
ERA=$2
FINAL_OUTDIR=$3

BASENAME=$(basename "$INPUT_FILE" .root)
JOB_OUTDIR="${FINAL_OUTDIR}/__job_${BASENAME}_${ERA}"
mkdir -p "$JOB_OUTDIR"

python3 run_processor.py --input "$INPUT_FILE" --output "$JOB_OUTDIR" --dataset "$ERA" --year 2024 --metaconditions /eos/home-a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/correctionsAndPresel/metaconditions.json --corrections "Scale_IJazZ" --treename Events

for pq in "$JOB_OUTDIR"/photons*.parquet; do
    [ -e "$pq" ] || continue
    VARIATION=$(basename "$pq" .parquet | sed "s/^photons//")
    CF="$JOB_OUTDIR/cutflow${VARIATION}.csv"
    mv "$pq" "${FINAL_OUTDIR}/${BASENAME}${VARIATION}_photons.parquet"
    if [ -e "$CF" ]; then
        mv "$CF" "${FINAL_OUTDIR}/${BASENAME}${VARIATION}_cutflow.csv"
    fi
done

rm -rf "$JOB_OUTDIR"
