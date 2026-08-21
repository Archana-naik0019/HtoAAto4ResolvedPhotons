#!/bin/bash
JOB_IDX=$1
JOB_DATA=$(python3 -c "import json, base64; data=json.load(open('/eos/home-a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/NanoAOD_skimming/Skimmer/condor_jobs_skimmer/job_map.json')); item=data[$JOB_IDX]; print(base64.b64decode(item['encoded_input_file']).decode('utf-8'), item['dataset'], item['outdir'], item['outfile'])")
read -r INPUT_FILE DATASET JOB_OUTDIR OUT_FILE <<< "$JOB_DATA"

mkdir -p "$JOB_OUTDIR"
exec > "$JOB_OUTDIR/job.out" 2> "$JOB_OUTDIR/job.err"

set -e

source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc15-opt/setup.sh
export PYTHONPATH="$HOME/.local/lib/python3.13/site-packages/:$PYTHONPATH"
export PATH="$HOME/.local/bin/:$PATH"

export SETUPTOOLS_SCM_PRETEND_VERSION="22.0.0"
export SETUPTOOLS_SCM_PRETEND_VERSION_FOR_PYARROW="22.0.0"

cd /eos/home-a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/NanoAOD_skimming/Skimmer

python3 /eos/home-a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/NanoAOD_skimming/Skimmer/nano_reduce.py --input "$INPUT_FILE" --output "$OUT_FILE" --is-data --lumimask "/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/files_for_evtMixing/Cert_Collisions2024_378981_386951_Golden.json" 
