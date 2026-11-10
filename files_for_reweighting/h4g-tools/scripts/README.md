# BDT & AMS Machinery

This directory contains scripts for BDT training, scoring, and AMS (category) calculation.

## Workflow Commands

### 1. BDT Training
Generate and train the BDT model:
```bash
python3 patch_run_fix.py --gen -fp -tj <JSON file name with input parameters in /h4g-tools/scripts/bdtIO/> -sig <dir with processed signal samples> -y 2024 -xgb <model name> > gen_bdt.log 2>&1
```

BDT Scoring/Inferencing:
```bash
python3 patch_run_fix.py --run -sig <dir with processed signal samples> -y 2024 -xgb <model name> > run_bdt.log 2>&1
```

Category Calculation (defines the BDT cuts based on AMS maximization):
```bash
python3 patch_cats.py -y 2024 -xgb <model name> --nCats 1 > calc_categories.log 2>&1
```

Applying the obtained BDT cuts:
```bash
python3 patch_applyCuts.py -y 2024 -xgb <model name> --nCats 1 > apply_categories.log 2>&1
```


