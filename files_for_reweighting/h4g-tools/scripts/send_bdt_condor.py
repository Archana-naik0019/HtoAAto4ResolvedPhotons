import os
import argparse
import subprocess
from time import strftime
from h4g_tools.utils.runner_utils import get_parser_scan, handle_errors

##  Get directory of file
cwd = os.path.dirname(os.path.abspath(__file__))


if __name__ == "__main__":
    # Command-line arguments setup
    parser = get_parser_scan()
    parser.add_argument("-q", "--queue", default="microcentury", choices=["espresso", "microcentury", "longlunch", "workday", "tomorrow", "testmatch", "nextweek"], required=False, help="Specify queue for resubmission. Defaults to microcentury.")
    parser.add_argument("-ram", "--ram", type=int, required=False, default=24, help="Specify RAM usage for submission. Defaults to 24GB. Look at ResidentSetSize for maximum RAM used and set the jobs sligtly above that number.")
    parser.add_argument("-tj", "--training_json", type=str, required=False, help="Provide a training json for input parameters. Relative path to bdtIO.")
    args = parser.parse_args()

    # Catch parser issues
    handle_errors(
        ["Need training samples.", args.gen, args.sigInput is None],
        ["Only need to provide signal samples and run reweighting before this.", args.gen, args.bkgInput is not None, args.data is not None],
        ["Need samples to run over.", args.run, args.data is None, args.sigInput is None, args.bkgInput is None],
        ["Cannot debug when running or plotting BDT results.", args.debug, args.plot or args.run],
        ["Must provide a year for processing.", args.year is None]
    )

    submission_dir = f"bdt_{args.year}_{args.xgb_name}_" + subprocess.getoutput("date +%Y%m%d_%H%M%S")

    submit_file = f"""executable = {cwd}/.bdt_submissions/{submission_dir}/execs/{{0}}
arguments = $(Process)
output = {{1}}/$(ClusterId).$(Process).out
error = {{1}}/$(ClusterId).$(Process).err
log = {{1}}/$(ClusterId).log
request_memory = {{2}}GB
batch_name = bdt{strftime("%d%H%M")}
getenv = True
+JobFlavour = "{{3}}"
on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)
periodic_remove = JobStatus == 1 && CurrentTime-EnteredCurrentStatus > 3600
max_retries = 2
queue {{4}}
"""

    sh_file = """#!/bin/sh
unset PYTHONPATH
# No VOMS proxy needed!
"""

    print(f"Writing jobs for BDT trained on {args.training_json}...")

    # Make job submission directories!
    if not os.path.exists(os.path.join(cwd, ".bdt_submissions")):
        os.mkdir(os.path.join(cwd, ".bdt_submissions"))
    if not os.path.exists(os.path.join(cwd, ".bdt_submissions", submission_dir)):
        os.mkdir(os.path.join(cwd, ".bdt_submissions", submission_dir))
    if not os.path.exists(os.path.join(cwd, ".bdt_submissions", submission_dir, "jobs")):
        os.mkdir(os.path.join(cwd, ".bdt_submissions", submission_dir, "jobs"))
    if not os.path.exists(os.path.join(cwd, ".bdt_submissions", submission_dir, "execs")):
        os.mkdir(os.path.join(cwd, ".bdt_submissions", submission_dir, "execs"))
    if not os.path.exists(os.path.join(cwd, ".bdt_submissions", submission_dir, "logs")):
        os.mkdir(os.path.join(cwd, ".bdt_submissions", submission_dir, "logs"))

    # Write one sh and sub file for all default weights. Each job will then get ONE queue number to run ONE of the if-statements in the sh file per file.
    sub_path = f"condor_bdt.sub"
    sh_path = f"exec_bdt.sh"

    # One sh file for all jobs (files) in a given era
    with open(os.path.join(cwd, ".bdt_submissions", submission_dir, "execs", sh_path), "w") as executable:
        executable.write(sh_file)
        executable.write("\n")

        # Start if-statementt
        executable.write(f"if [ $1 -eq 0 ]; then\n")

        # Write bdt command
        executable.write(f" /usr/bin/env python3 {os.path.abspath(os.path.join(cwd, 'scan_masses_bdt.py'))} --gen -fp -sig {args.sigInput} -y {args.year} -xgb {args.xgb_name} -tj {args.training_json}".strip())
        executable.write("\n")
        
        # Finish if-statement
        executable.write(f"fi\n\n")

    with open(os.path.join(cwd, ".bdt_submissions", submission_dir, "jobs", sub_path), "w") as submission:
        submission.write(submit_file.format(
            sh_path,
            os.path.join(cwd, ".bdt_submissions", submission_dir, "logs"),
            args.ram,
            args.queue,
            1
        ))

    # Set permissions of executable so condor can use it
    assert os.path.exists(os.path.join(cwd, ".bdt_submissions", submission_dir, "execs", sh_path))
    os.system(f"chmod 775 {os.path.join(cwd, '.bdt_submissions', submission_dir, 'execs', sh_path)}")

    # Submit jobs
    assert os.path.exists(os.path.join(cwd, ".bdt_submissions", submission_dir, "jobs", sub_path))
    subprocess.run(["condor_submit", os.path.join(cwd, ".bdt_submissions", submission_dir, "jobs", sub_path)])


