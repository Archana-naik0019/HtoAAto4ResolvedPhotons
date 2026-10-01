# H4g Tools
This repository is for tools for the H4g analysis including:
 * BDT training and running

```bash
# H4g Tools Setup
git clone ssh://git@gitlab.cern.ch:7999/castells/h4g-tools.git
cd h4g-tools/

conda create --file environment.yml  # I would recommend mamba for this...
conda activate h4g-tools  # make empty environment if necessary to contain everything nicely
python3 -m pip install -e .
```

NOTE:
Make sure you clone paper_plots into the scripts directory. The plotting script (generate_plots.py) dumps a ROOT files with histograms there.
