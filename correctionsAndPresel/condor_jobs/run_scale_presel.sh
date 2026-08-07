#!/bin/bash

source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc15-opt/setup.sh
cd /eos/home-a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/correctionsAndPresel

python3 scale_presel_bulk_submissions.py --input-dir "$1" --output-dir "$2" --metaconditions metaconditions.json --year 2024 --corrections "Scale_IJazZ" --treename Events
