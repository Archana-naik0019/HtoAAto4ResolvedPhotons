#!/usr/bin/env python3

import subprocess
import json
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input", help="Text file containing DAS dataset names")
parser.add_argument("output", help="Output JSON")
parser.add_argument(
    "--redirector",
    default="root://cms-xrd-global.cern.ch/",
    help="XRootD redirector"
)
args = parser.parse_args()

sample_dict = {}

with open(args.input) as f:
    datasets = [line.strip() for line in f if line.strip()]

for dataset in datasets:

    print(f"Processing {dataset}")

    # Extract era (Run2024C, Run2024D, ...)
    m = re.search(r"Run2024([A-Z])", dataset)
    if m:
        era = "2024" + m.group(1)
    else:
        era = "Unknown"

    if era not in sample_dict:
        sample_dict[era] = []

    cmd = f'dasgoclient -query="file dataset={dataset}"'
    result = subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print(f"Failed: {dataset}")
        continue

    files = result.stdout.strip().split("\n")

    for file in files:
        if file:
            sample_dict[era].append(args.redirector + file)

with open(args.output, "w") as f:
    json.dump(sample_dict, f, indent=4)

print(f"\nWritten {args.output}")
