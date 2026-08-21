#!/usr/bin/env python3
"""
condorSubmit_skimmer.py

Submits HTCondor jobs per file using base64 URL encoding to avoid spool parser crashes,
redirecting stdout/stderr directly to EOS inside the bash wrapper to execute nano_reduce.py.
"""

import argparse
import base64
import json
import os
import pathlib
import subprocess


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "--fileset",
        default="/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/NanoAOD_skimming/Skimmer/2024C_6failed.json",
        help="Path to input fileset JSON: {dataset_name: [file_urls]}",
    )
    parser.add_argument(
        "--output-dir",
        default="/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/NanoAOD_skimming/Skimmer/2024C_updated6missing",
        help="Top-level output dir on EOS",
    )
    parser.add_argument(
        "--lumimask",
        default="/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/files_for_evtMixing/Cert_Collisions2024_378981_386951_Golden.json",
        help="Path to Golden JSON file for lumi mask filtering",
    )
    parser.add_argument(
        "--is-data",
        action="store_true",
        default=True,
        help="Flag to treat inputs as Data (enables Lumi Mask and Trigger filtering)",
    )
    parser.add_argument(
        "--reducer-script",
        default="nano_reduce.py",
        help="Path to the nano_reduce.py script",
    )

    # ---- HTCondor knobs ----
    parser.add_argument(
        "--jobflavour",
        default="nextweek",
        help="espresso, microcentury, longlunch, workday, tomorrow, testmatch, nextweek",
    )
    parser.add_argument(
        "--memory", default="30000MB", help="Per-job memory request"
    )
    parser.add_argument(
        "--lcg-view",
        default="/cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc15-opt/setup.sh",
        help="LCG environment setup script sourced by the job wrapper",
    )

    return parser


def format_eos_root_url(path):
    """Converts local /eos/... paths to root://eosuser.cern.ch//eos/... URL format for CERN Condor."""
    abs_path = os.path.abspath(path)
    if abs_path.startswith("/eos/"):
        return f"root://eosuser.cern.ch/{abs_path}"
    return abs_path


def get_x509_proxy():
    """Finds active x509 user proxy dynamically."""
    if "X509_USER_PROXY" in os.environ and os.path.exists(
        os.environ["X509_USER_PROXY"]
    ):
        return os.environ["X509_USER_PROXY"]

    uid = os.getuid()
    default_proxy = f"/tmp/x509up_u{uid}"
    if os.path.exists(default_proxy):
        return default_proxy

    return None


def main():
    args = get_parser().parse_args()

    submit_dir = "condor_jobs_skimmer"
    os.makedirs(submit_dir, exist_ok=True)
    os.makedirs(os.path.join(submit_dir, "logs"), exist_ok=True)

    workdir_abs = os.getcwd()
    lumimask_abs = os.path.abspath(args.lumimask) if args.lumimask else None
    reducer_script_abs = os.path.abspath(args.reducer_script)

    proxy_path = get_x509_proxy()

    with open(args.fileset) as f:
        fileset = json.load(f)

    output_dir_abs = os.path.abspath(args.output_dir)
    os.makedirs(output_dir_abs, exist_ok=True)

    job_map = []
    for dataset_name, file_list in fileset.items():
        for file_url in file_list:
            file_tag = pathlib.Path(file_url).stem
            job_outdir = os.path.join(output_dir_abs, dataset_name, file_tag)
            os.makedirs(job_outdir, exist_ok=True)

            out_root_file = os.path.join(job_outdir, f"skim_{file_tag}.root")

            encoded_url = base64.b64encode(file_url.encode("utf-8")).decode(
                "utf-8"
            )

            job_map.append(
                {
                    "encoded_input_file": encoded_url,
                    "dataset": dataset_name,
                    "outdir": job_outdir,
                    "outfile": out_root_file,
                }
            )

    job_map_file = os.path.abspath(os.path.join(submit_dir, "job_map.json"))
    with open(job_map_file, "w") as f:
        json.dump(job_map, f, indent=2)

    njobs = len(job_map)

    ###############################################################
    # Shell Wrapper: Executes nano_reduce.py per file
    ###############################################################
    shell = os.path.join(submit_dir, "run_single_file.sh")

    with open(shell, "w") as f:
        f.write("#!/bin/bash\n")

        f.write("JOB_IDX=$1\n")
        f.write(
            f'JOB_DATA=$(python3 -c "import json, base64; data=json.load(open(\'{job_map_file}\')); item=data[$JOB_IDX]; print(base64.b64decode(item[\'encoded_input_file\']).decode(\'utf-8\'), item[\'dataset\'], item[\'outdir\'], item[\'outfile\'])")\n'
        )
        f.write('read -r INPUT_FILE DATASET JOB_OUTDIR OUT_FILE <<< "$JOB_DATA"\n\n')

        f.write('mkdir -p "$JOB_OUTDIR"\n')
        f.write('exec > "$JOB_OUTDIR/job.out" 2> "$JOB_OUTDIR/job.err"\n\n')
        f.write("set -e\n\n")

        # 1. Source LCG stack
        f.write(f"source {args.lcg_view}\n")

        # 2. Append local user packages AFTER CVMFS is sourced
        f.write('export PYTHONPATH="$HOME/.local/lib/python3.13/site-packages/:$PYTHONPATH"\n')
        f.write('export PATH="$HOME/.local/bin/:$PATH"\n\n')

        # 3. Setup setuptools_scm pretend versions
        f.write('export SETUPTOOLS_SCM_PRETEND_VERSION="22.0.0"\n')
        f.write('export SETUPTOOLS_SCM_PRETEND_VERSION_FOR_PYARROW="22.0.0"\n\n')

        f.write(f"cd {workdir_abs}\n\n")

        # Build nano_reduce.py command
        cmd = f'python3 {reducer_script_abs} --input "$INPUT_FILE" --output "$OUT_FILE" '
        if args.is_data:
            cmd += "--is-data "
        if lumimask_abs:
            cmd += f'--lumimask "{lumimask_abs}" '
        cmd += "\n"

        f.write(cmd)

    os.chmod(shell, 0o755)

    ###############################################################
    # HTCondor Job Index List
    ###############################################################
    jobs = os.path.join(submit_dir, "jobs.txt")
    with open(jobs, "w") as f:
        for idx in range(njobs):
            f.write(f"{idx}\n")

    ###############################################################
    # HTCondor Submit Description File
    ###############################################################
    submit = os.path.join(submit_dir, "submit.sub")

    with open(submit, "w") as f:
        f.write(f"""universe = vanilla

executable = {os.path.abspath(shell)}
arguments = $(job_idx)

""")
        if proxy_path:
            f.write(f"x509userproxy = {proxy_path}\n")
            f.write("use_x509userproxy = True\n")
            f.write(f"transfer_input_files = {proxy_path}\n")

        f.write(f"""
notification = never
should_transfer_files = YES
when_to_transfer_output = ON_EXIT_OR_EVICT
MY.OutputDestination = "{format_eos_root_url(output_dir_abs)}"
getenv = True

request_memory = {args.memory}
request_cpus = 1

+JobFlavour = "{args.jobflavour}"

log = {os.path.abspath(submit_dir)}/logs/$(Cluster).log

queue job_idx from {os.path.abspath(jobs)}
""")

    print(
        f"\nCreated {njobs} job(s), one per ROOT file, across {len(fileset)} dataset(s):"
    )
    for dataset_name, file_list in fileset.items():
        print(f"  {dataset_name}: {len(file_list)} file(s)")

    print("\nSubmitting jobs to HTCondor...")
    subprocess.run(["condor_submit", "--spool", submit], check=True)


if __name__ == "__main__":
    main()
