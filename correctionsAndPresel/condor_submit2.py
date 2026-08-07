#!/usr/bin/env python3
"""
condor_submitter.py

Submits one HTCondor job per ROOT file, each directly invoking
run_processor.py (no scale_presel_bulk_submissions.py wrapper).

Each job:
  - runs in a private scratch output dir
  - run_processor.py writes photons*.parquet / cutflow*.csv there
  - the wrapper script renames/moves them into
      <output-dir>/<era>/<basename>_photons.parquet
      <output-dir>/<era>/<basename>_cutflow.csv
    (handles either "photons.parquet" or "photons_<variation>.parquet"
    naming, whichever your run_processor.py currently produces)
"""

import argparse
import glob
import os
import subprocess


def get_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input-dir", required=True,
                         help="Top-level dir containing era folders (2024C, 2024D, ...)")
    parser.add_argument("--output-dir", required=True,
                         help="Top-level output dir (EOS is fine here)")
    parser.add_argument("--metaconditions", required=True)

    parser.add_argument("--year", default="2024")
    parser.add_argument("--corrections", default="")
    parser.add_argument("--systematics", default="")
    parser.add_argument("--treename", default="Events")

    parser.add_argument("--doFlow_corrections", action="store_true")
    parser.add_argument("--use_photonid_mva", action="store_true")

    parser.add_argument(
        "--jobflavour",
        default="nextweek",
        help="espresso, microcentury, longlunch, workday, tomorrow, testmatch, nextweek",
    )
    parser.add_argument("--memory", default="4000MB")

    return parser


def main():
    args = get_parser().parse_args()

    submit_dir = "condor_jobs"
    os.makedirs(submit_dir, exist_ok=True)
    os.makedirs(os.path.join(submit_dir, "logs"), exist_ok=True)
    os.makedirs(os.path.join(submit_dir, "out"), exist_ok=True)
    os.makedirs(os.path.join(submit_dir, "err"), exist_ok=True)

    metaconditions_abs = os.path.abspath(args.metaconditions)
    workdir_abs = os.getcwd()

    ###############################################################
    # shell wrapper - runs run_processor.py directly, then renames
    # its output into the flat <era>/<basename>_photons.parquet scheme
    ###############################################################
    shell = os.path.join(submit_dir, "run_single_file.sh")

    with open(shell, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("set -e\n\n")
        f.write("source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc15-opt/setup.sh\n")
        f.write(f"cd {workdir_abs}\n\n")

        f.write("INPUT_FILE=$1\n")
        f.write("ERA=$2\n")
        f.write("FINAL_OUTDIR=$3\n\n")

        f.write('BASENAME=$(basename "$INPUT_FILE" .root)\n')
        f.write('JOB_OUTDIR="${FINAL_OUTDIR}/__job_${BASENAME}_${ERA}"\n')
        f.write('mkdir -p "$JOB_OUTDIR"\n\n')

        cmd = (
            'python3 run_processor.py '
            '--input "$INPUT_FILE" '
            '--output "$JOB_OUTDIR" '
            '--dataset "$ERA" '
            f'--year {args.year} '
            f'--metaconditions {metaconditions_abs} '
        )
        if args.corrections:
            cmd += f'--corrections "{args.corrections}" '
        if args.systematics:
            cmd += f'--systematics "{args.systematics}" '
        if args.doFlow_corrections:
            cmd += "--doFlow_corrections "
        if args.use_photonid_mva:
            cmd += "--use_photonid_mva "
        cmd += f"--treename {args.treename}\n\n"
        f.write(cmd)

        # Move/rename whatever photons*.parquet + matching cutflow*.csv
        # were produced, into the flat final naming, then clean up.
        f.write('for pq in "$JOB_OUTDIR"/photons*.parquet; do\n')
        f.write('    [ -e "$pq" ] || continue\n')
        f.write('    VARIATION=$(basename "$pq" .parquet | sed "s/^photons//")\n')
        f.write('    CF="$JOB_OUTDIR/cutflow${VARIATION}.csv"\n')
        f.write('    mv "$pq" "${FINAL_OUTDIR}/${BASENAME}${VARIATION}_photons.parquet"\n')
        f.write('    if [ -e "$CF" ]; then\n')
        f.write('        mv "$CF" "${FINAL_OUTDIR}/${BASENAME}${VARIATION}_cutflow.csv"\n')
        f.write('    fi\n')
        f.write('done\n\n')
        f.write('rm -rf "$JOB_OUTDIR"\n')

    os.chmod(shell, 0o755)

    ###############################################################
    # jobs.txt - one line per ROOT file
    ###############################################################
    jobs = os.path.join(submit_dir, "jobs.txt")
    njobs = 0

    with open(jobs, "w") as f:
        eras = sorted(
            d for d in os.listdir(args.input_dir)
            if os.path.isdir(os.path.join(args.input_dir, d))
        )
        for era in eras:
            files = sorted(glob.glob(os.path.join(args.input_dir, era, "*.root")))
            outdir = os.path.join(os.path.abspath(args.output_dir), era)
            os.makedirs(outdir, exist_ok=True)
            for infile in files:
                f.write(f"{infile} {era} {outdir}\n")
                njobs += 1

    ###############################################################
    # submit file
    ###############################################################
    submit = os.path.join(submit_dir, "submit.sub")

    with open(submit, "w") as f:
        f.write(f"""universe = vanilla

executable = {os.path.abspath(shell)}
arguments = $(input) $(era) $(outdir)

notification = never
should_transfer_files = YES
when_to_transfer_output = ON_EXIT_OR_EVICT
transfer_output_files = ""
getenv = True

request_memory = {args.memory}
request_cpus = 1

+JobFlavour = "{args.jobflavour}"

output = {submit_dir}/out/$(Cluster)_$(Process).out
error  = {submit_dir}/err/$(Cluster)_$(Process).err
log    = {submit_dir}/logs/$(Cluster).log

queue input,era,outdir from {os.path.abspath(jobs)}
""")

    print(f"\nCreated {njobs} jobs")

    subprocess.run(
        ["condor_submit", "--spool", submit],
        check=True,
    )


if __name__ == "__main__":
    main()
