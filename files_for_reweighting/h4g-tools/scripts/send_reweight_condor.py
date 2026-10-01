import os
import argparse
import subprocess
from time import strftime
from h4g_tools.utils.runner_utils import handle_errors

##  Get directory of file
cwd = os.path.dirname(os.path.abspath(__file__))


if __name__ == "__main__":
    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser", conflict_handler="resolve")
    parser.add_argument("-bkg", "--bkgInput", help="Directory to import background samples from.", type=str, default=None, required=True)
    parser.add_argument("-d", "--data", help="Directory to import data samples from.", type=str, default=None, required=True)
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    parser.add_argument("-nb", "--nBins", type=int, required=False, help="Number of bins to use (excluding under/overflow).")
    parser.add_argument("-q", "--queue", default="longlunch", choices=["espresso", "microcentury", "longlunch", "workday", "tomorrow", "testmatch", "nextweek"], required=False, help="Specify queue for resubmission. Defaults to microcentury.")
    parser.add_argument("-ram", "--ram", type=int, required=False, default=8, help="Specify RAM usage for submission. Defaults to 8GB. Look at ResidentSetSize for maximum RAM used and set the jobs sligtly above that number.")

    args = parser.parse_args()
    
    # Catch parser issues
    handle_errors(
        ["Need samples to run over.", args.data is None and args.bkgInput is None],
        ["Must provide a year for processing.", args.year is None],
    )

    submission_dir = "reweighting_" + args.year + "_" + subprocess.getoutput("date +%Y%m%d_%H%M%S")

    submit_file = f"""executable = {cwd}/.reweighting_submissions/{submission_dir}/execs/{{0}}
arguments = $(Process)
output = {{1}}/$(ClusterId).$(Process).out
error = {{1}}/$(ClusterId).$(Process).err
log = {{1}}/$(ClusterId).log
request_memory = {{2}}GB
batch_name = reweighting{strftime("%d%H%M")}
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

    #default_weights = ["1.0", "0.75", "0.50", "0.25", "0.1", "0.05", "0.015", "0.025"]
    default_weights = ["0.1"]
    print("Writing jobs for default weights: ", ", ".join([str(dw) for dw in default_weights]))

    # Make job submission directories!
    if not os.path.exists(os.path.join(cwd, ".reweighting_submissions")):
        os.mkdir(os.path.join(cwd, ".reweighting_submissions"))
    if not os.path.exists(os.path.join(cwd, ".reweighting_submissions", submission_dir)):
        os.mkdir(os.path.join(cwd, ".reweighting_submissions", submission_dir))
    if not os.path.exists(os.path.join(cwd, ".reweighting_submissions", submission_dir, "jobs")):
        os.mkdir(os.path.join(cwd, ".reweighting_submissions", submission_dir, "jobs"))
    if not os.path.exists(os.path.join(cwd, ".reweighting_submissions", submission_dir, "execs")):
        os.mkdir(os.path.join(cwd, ".reweighting_submissions", submission_dir, "execs"))
    if not os.path.exists(os.path.join(cwd, ".reweighting_submissions", submission_dir, "logs")):
        os.mkdir(os.path.join(cwd, ".reweighting_submissions", submission_dir, "logs"))

    # Write one sh and sub file for all default weights. Each job will then get ONE queue number to run ONE of the if-statements in the sh file per file.
    sub_path = f"condor_reweight.sub"
    sh_path = f"exec_reweight.sh"

    # One sh file for all jobs (files) in a given era
    with open(os.path.join(cwd, ".reweighting_submissions", submission_dir, "execs", sh_path), "w") as executable:
        executable.write(sh_file)
        executable.write("\n")

        for idx, dw in enumerate(default_weights):
            # Start if-statementt
            executable.write(f"if [ $1 -eq {idx} ]; then\n")

            # Write reweighting command
            executable.write((f" /usr/bin/env python3 {os.path.abspath(os.path.join(cwd, 'compute_Ndim_reweighting.py'))} -bkg {args.bkgInput} -d {args.data} -dw -dwf {dw} -y {args.year}" + ((f" -nb {args.nBins}") if args.nBins is not None else "")).strip())
            executable.write("\n")
            
            # Finish if-statement
            executable.write(f"fi\n\n")

    with open(os.path.join(cwd, ".reweighting_submissions", submission_dir, "jobs", sub_path), "w") as submission:
        submission.write(submit_file.format(
            sh_path,
            os.path.join(cwd, ".reweighting_submissions", submission_dir, "logs"),
            args.ram,
            args.queue,
            len(default_weights)
        ))

    # Set permissions of executable so condor can use it
    assert os.path.exists(os.path.join(cwd, ".reweighting_submissions", submission_dir, "execs", sh_path))
    os.system(f"chmod 775 {os.path.join(cwd, '.reweighting_submissions', submission_dir, 'execs', sh_path)}")

    # Submit jobs
    assert os.path.exists(os.path.join(cwd, ".reweighting_submissions", submission_dir, "jobs", sub_path))
    subprocess.run(["condor_submit", os.path.join(cwd, ".reweighting_submissions", submission_dir, "jobs", sub_path)])


