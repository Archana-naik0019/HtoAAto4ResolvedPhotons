import os
import argparse
import subprocess
from time import strftime
from h4g_tools.utils.runner_utils import get_parser_cats, handle_errors

##  Get directory of file
cwd = os.path.dirname(os.path.abspath(__file__))


if __name__ == "__main__":
    # Command-line arguments setup
    parser = get_parser_cats()
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    parser.add_argument("-xgb", "--xgb_name", type=str, required=True, help="XGBoost model name.")
    parser.add_argument("-q", "--queue", default="microcentury", choices=["espresso", "microcentury", "longlunch", "workday", "tomorrow", "testmatch", "nextweek"], required=False, help="Specify queue for resubmission. Defaults to microcentury.")
    parser.add_argument("-ram", "--ram", type=int, required=False, default=24, help="Specify RAM usage for submission. Defaults to 24GB. Look at ResidentSetSize for maximum RAM used and set the jobs sligtly above that number.")
    args = parser.parse_args()

    handle_errors(
        ["Span must be between 0.0 and 1.0", (args.span < 0.0 or args.span > 1.0), (args.span != -1.0)]
    )

    submission_dir = f"categories_{args.year}_{args.nCats}_" + subprocess.getoutput("date +%Y%m%d_%H%M%S")

    submit_file = f"""executable = {cwd}/.categories_submissions/{submission_dir}/execs/{{0}}
arguments = $(Process)
output = {{1}}/$(ClusterId).$(Process).out
error = {{1}}/$(ClusterId).$(Process).err
log = {{1}}/$(ClusterId).log
request_memory = {{2}}GB
batch_name = categories{strftime("%d%H%M")}
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

    print(f"Writing jobs for {args.nCats} categories...")

    # Make job submission directories!
    if not os.path.exists(os.path.join(cwd, ".categories_submissions")):
        os.mkdir(os.path.join(cwd, ".categories_submissions"))
    if not os.path.exists(os.path.join(cwd, ".categories_submissions", submission_dir)):
        os.mkdir(os.path.join(cwd, ".categories_submissions", submission_dir))
    if not os.path.exists(os.path.join(cwd, ".categories_submissions", submission_dir, "jobs")):
        os.mkdir(os.path.join(cwd, ".categories_submissions", submission_dir, "jobs"))
    if not os.path.exists(os.path.join(cwd, ".categories_submissions", submission_dir, "execs")):
        os.mkdir(os.path.join(cwd, ".categories_submissions", submission_dir, "execs"))
    if not os.path.exists(os.path.join(cwd, ".categories_submissions", submission_dir, "logs")):
        os.mkdir(os.path.join(cwd, ".categories_submissions", submission_dir, "logs"))

    # Write one sh and sub file for all default weights. Each job will then get ONE queue number to run ONE of the if-statements in the sh file per file.
    sub_path = f"condor_categories.sub"
    sh_path = f"exec_categories.sh"

    # One sh file for all jobs (files) in a given era
    with open(os.path.join(cwd, ".categories_submissions", submission_dir, "execs", sh_path), "w") as executable:
        executable.write(sh_file)
        executable.write("\n")

        # Start if-statementt
        executable.write(f"if [ $1 -eq 0 ]; then\n")

        # Write categories command
        executable.write(f" /usr/bin/env python3 {os.path.abspath(os.path.join(cwd, 'calculate_categories.py'))} -y {args.year} -xgb {args.xgb_name} -c {args.nCats}".strip())
        executable.write("\n")
        
        # Finish if-statement
        executable.write(f"fi\n\n")

    with open(os.path.join(cwd, ".categories_submissions", submission_dir, "jobs", sub_path), "w") as submission:
        submission.write(submit_file.format(
            sh_path,
            os.path.join(cwd, ".categories_submissions", submission_dir, "logs"),
            args.ram,
            args.queue,
            1
        ))

    # Set permissions of executable so condor can use it
    assert os.path.exists(os.path.join(cwd, ".categories_submissions", submission_dir, "execs", sh_path))
    os.system(f"chmod 775 {os.path.join(cwd, '.categories_submissions', submission_dir, 'execs', sh_path)}")

    # Submit jobs
    assert os.path.exists(os.path.join(cwd, ".categories_submissions", submission_dir, "jobs", sub_path))
    subprocess.run(["condor_submit", os.path.join(cwd, ".categories_submissions", submission_dir, "jobs", sub_path)])


