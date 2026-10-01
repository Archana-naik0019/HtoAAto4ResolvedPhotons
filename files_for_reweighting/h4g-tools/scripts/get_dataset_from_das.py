import os
from pathlib import Path
import json
import argparse
import subprocess
import shlex

"""
Example:
python3 get_dataset_from_das.py datasets.txt
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Command line options parser")
    parser.add_argument("dataset_file", help="Dataset file to load.", type=str)
    parser.add_argument("-of", "--outputFile", help="Output JSON for datasets.", type=str, default=None, required=False)
    args = parser.parse_args()

    ## Output file setup
    if args.outputFile is None:
        args.outputFile = "tmp_dataset.json"

    print(f"Using output filename: {args.outputFile}")
    if os.path.exists(args.outputFile):
        assert not os.path.exists(args.outputFile), "Output file will not be overwritten! Backup the file and try again."

    ## Get datasets from file
    dbs_out = []
    datasets = []
    with open(args.dataset_file, "r") as f:
        for line in f:
            datasets.append(line.strip())

    ## Get files in datasets from DAS using dasgoclient
    for dataset in datasets:
        dbs_out.extend(
            subprocess.run(
                shlex.split(f'dasgoclient -query "file dataset={dataset}"'),
                stdout = subprocess.PIPE,
                universal_newlines = True).stdout.split("\n")[:-1]
        )

    ## Fill dictionary with file paths
    d = {}
    previous_egamma = ""
    current_egamma = ""
    previous_dataset = ""
    current_dataset = ""
    # EDIT THE REDICTOR AS YOU NEED!!!!
    redirector = "root://cmsxcache.crc.nd.edu//"
    for line in dbs_out:
        # Change this to 4 for H4g signal MC. For EGamma data this works with 3.
        era = None
        offset = 3
        if "Run3Summer22" in line:
            if "postEE" in line:
                era = "_postEE"
            else:
                era = "_preEE"
        elif "Run3Summer23" in line:
            if "BPix" in line:
                era = "_postBPix"
            else:
                era = "_preBPix"
        elif "RunIII2024" in line:
            era = ""

        if era is not None:
            # For signal MC and bkg MC
            offset = 4

        egamma_offset = offset
        if "2022_dataset" in args.dataset_file:
            egamma_offset -= 1
        current_egamma = line.split("/")[egamma_offset+1]
        current_dataset = line.split("/")[offset]
        if era is None:
            if current_egamma.lower() != previous_egamma.lower():
                print(f"Handling {current_egamma} for {current_dataset}...")
        
        if era is not None and ("GG" not in line and "GJet" not in line and "QCD" not in line):
            current_dataset = "Signal_" + current_dataset[current_dataset.find("MA-")+3:current_dataset.find("MA-")+5] + "_GeV" + era
        elif "GG" in line or "GJet" in line or "QCD" in line:
            pass

        if current_dataset.lower() != previous_dataset.lower():
            d.update({current_dataset: list().copy()})

        #print(offset, era, current_dataset, line)
        #break

        d[current_dataset].append(f"{redirector}{line}")
        previous_egamma = current_egamma
        previous_dataset = current_dataset

    ## Dump datasets/files dict to JSON
    with open(args.outputFile, "w") as f:
        json.dump(
            d,
            f,
            indent = 4,
            sort_keys=True
        )
        print(f"Dumped JSON to {args.outputFile}!")
