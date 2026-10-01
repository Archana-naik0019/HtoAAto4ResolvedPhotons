import os
from pathlib import Path
import json
import argparse
import subprocess
import shlex

"""
Make sure to have a GRID proxy and to load brilcalc first:
voms-proxy-init --voms cms
source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env
export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda3/bin:$PATH

Example:
lumi_from_das.py -d /EGamma/Run2018A-UL2018_MiniAODv2-v1/MINIAOD,/EGamma/Run2018B-UL2018_MiniAODv2-v1/MINIAOD,/EGamma/Run2018C-UL2018_MiniAODv2-v1/MINIAOD,/EGamma/Run2018D-UL2018_MiniAODv2-v2/MINIAOD -of test.txt
"""

class splitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        setattr(
            namespace,
            self.dest,
            values.split(",")
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Command line options parser")
    parser.add_argument("-d", "--dataset", help="Dataset(s) to query", action=splitArgs, type=str, required=True)
    parser.add_argument("-of", "--outputFile", help="Output file", type=str, default=None, required=False)
    parser.add_argument("-hlt", "--hlt", help="HLT", type=str, default="", required=False)
    args = parser.parse_args()

    if args.outputFile is None:
        args.outputFile = "tmp_lumi.txt"
    if os.path.exists(args.outputFile):
        os.remove(args.outputFile)

    dbs_out = []
    for dataset in args.dataset:
        dbs_out.extend(
            subprocess.run(
                shlex.split(f'dasgoclient -query "run lumi dataset={dataset}"'),
                stdout = subprocess.PIPE,
                universal_newlines = True).stdout.split("\n")[:-1]
        )

    d = {}
    for line in dbs_out:
        l = line.split(" ")
        lumi = sorted([int(x) for x in l[1].replace("\n","").strip("[").strip("]").split(",")])

        new = []
        mx = lumi[0]
        curr = []
        for idx in range(len(lumi)):
            curr = [mx, lumi[idx]]
            if idx == len(lumi)-1:
                if len(lumi) > 1:
                    new.append(curr)
                if len(lumi) == 1:
                    new.append([lumi[0], lumi[0]])
                continue
            elif lumi[idx+1] - lumi[idx] > 1:
                new.append(curr)
                mx = lumi[idx+1]

        d.update({l[0]: new})

    with open(args.outputFile, "w") as f:
        json.dump(d, f)
        print("Dumped JSON")

    # Got normtag from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiLUM#How_to_calculate_the_luminosity
    #"/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json",
    #"/afs/crc.nd.edu/user/s/scastel2/Private/HiggsDNA-newest/higgs_dna/metaconditions/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json",
    print(
        subprocess.run([
            "brilcalc",
            "lumi",
            "--normtag",
            "/afs/crc.nd.edu/user/s/scastel2/Private/HiggsDNA-newest/higgs_dna/metaconditions/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json",
            "-u",
            "/fb",
            "-i",
            args.outputFile,
            "--hltpath",
            args.hlt
        ],
        stdout = subprocess.PIPE,
        universal_newlines = True).stdout
    )
