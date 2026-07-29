#!/bin/bash

# Usage:
# ./LumiCalculator.sh -y <year> -i <input_json> -o <output.csv> [-b <begin_run>] [-e <end_run>]

# Default values
DATATAG="online"
UNIT="/fb"
BEGIN_RUN=""
END_RUN=""
OUTPUT=""

# Parse arguments
while getopts "y:i:o:b:e:" opt; do
  case $opt in
    y) YEAR="$OPTARG" ;;
    i) JSON="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
    b) BEGIN_VAL="$OPTARG" ;;
    e) END_VAL="$OPTARG" ;;
    *) echo "Usage: $0 -y <year: 2022|2023|2024> -i <input_json> -o <output.csv> [-b <begin_run>] [-e <end_run>]"; exit 1 ;;
  esac
done

# Basic checks
if [[ -z "$YEAR" || -z "$JSON" || -z "$OUTPUT" ]]; then
  echo "Error: Year, input JSON file, and output CSV file are required."
  echo "Usage: $0 -y <year: 2022|2023|2024> -i <input_json> -o <output.csv> [-b <begin_run>] [-e <end_run>]"
  exit 1
fi

echo "DEBUG: YEAR=$YEAR"
echo "DEBUG: JSON=$JSON"
echo "DEBUG: OUTPUT=$OUTPUT"
echo "DEBUG: BEGIN_RUN=$BEGIN_VAL"
echo "DEBUG: END_RUN=$END_VAL"

# Activate brilcalc environment
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env
echo "BRIL Work Suite should now be available."

# Set normtag and datatag
case "$YEAR" in
  2022|2023)
    NORMTAG="/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json"
    DATATAG_FLAG=""
    ;;
  2024)
    NORMTAG="/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json"
    DATATAG_FLAG="--datatag $DATATAG"
    ;;
  *)
    echo "Error: Unsupported year '$YEAR'. Choose from 2022, 2023, or 2024."
    exit 1
    ;;
esac

# # Run range for 2023 only
# [[ "$YEAR" == "2023" && -n "$BEGIN_VAL" ]] && BEGIN_RUN="--begin $BEGIN_VAL"
# [[ "$YEAR" == "2023" && -n "$END_VAL" ]] && END_RUN="--end $END_VAL"

# Full brilcalc command using Singularity (bypassing alias)
BRILCALC_CMD="singularity -s exec --env PYTHONPATH=/home/bril/.local/lib/python3.10/site-packages /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-cloud/brilws-docker:latest brilcalc"

# Full lumi calculation command
CMD="$BRILCALC_CMD lumi --begin $BEGIN_VAL --end $END_VAL -i $JSON -u $UNIT --normtag $NORMTAG $DATATAG_FLAG --output-style csv"

echo "Running: $CMD"
eval "$CMD" > "$OUTPUT"

echo "Saved output to $OUTPUT"
