import argparse
import json
import os
import subprocess
from time import strftime
from pathlib import Path



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=str, help="Sample json to process with event mixing procedure.")
    parser.add_argument("dump", type=str, help="Output directory where to dump output files.")
    parser.add_argument("-g", "--goldenJSON", type=str, required=False, help="Golden JSON path. Only for data.")
    parser.add_argument("-xrd", "--xrootd", required=False, choices=["cms", "fermi", "infn", "global"], help="Choose different xrootd redirector for submission.")
    parser.add_argument("-q", "--queue", default="microcentury", choices=["espresso", "microcentury", "longlunch", "workday", "tomorrow", "testmatch", "nextweek"], required=False, help="Specify queue for resubmission. Defaults to microcentury.")
    parser.add_argument("-ram", "--ram", type=int, required=False, default=5, help="Specify RAM usage for submission. Defaults to 5GB. Look at ResidentSetSize for maximum RAM used and set the jobs sligtly above that number.")
    parser.add_argument("-nc", "--Ncycles", type=int, default=1, required=True, help="Event mixing cycles to run. Defaults to 1.")
    parser.add_argument("-offset", "--cycles_offset", type=int, default=0, required=False, help="Offset to use when running cycles for event mixing. Defaults to 0.")
    parser.add_argument("--do-3offset", action="store_true", help="Does offsets for cycles in sets of 3 instead of 1.")
    parser.add_argument("--no-event-mixing", action="store_true", help="Turns off event mixing procedure and only skims files. Useful for debugging or skimming data samples.")
    parser.add_argument("--debug", action="store_true", help="Run in debug mode. Limits to 20 events!")
    parser.add_argument("--show-all", action="store_true", help="Print all per-event information.")
    parser.add_argument("--noSkim", action="store_true", help="Do not skim samples prior to event mixing. Skims afterwards to reduce storage requirements.")
    # --show-all should be OFF when submitting Condor jobs since output files take up a ton of space.
    args = parser.parse_args()

    # Get path of this file
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Submission directory with timestamp to discern submissions with same input JSON
    submission_dir = args.file.replace(".json", "").strip("/").split("/")[-1] + "_" + subprocess.getoutput("date +%Y%m%d_%H%M%S")

    # Add Condor submission stuff here too
    ##  Write submit and sh files
    submit_file = f"""executable = /afs/cern.ch/user/a/arnaik/.event_mixing_submission/{submission_dir}/execs/{{0}}
arguments = $(Process)
output = {{2}}/{{1}}.$(ClusterId).$(Process).out
error = {{2}}/{{1}}.$(ClusterId).$(Process).err
log = {{2}}/{{1}}.$(ClusterId).log
request_memory = {{3}}GB
batch_name = evtmix{strftime("%d%H%M")}
getenv = True
+JobFlavour = "{{4}}"
on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)
periodic_remove = JobStatus == 1 && CurrentTime-EnteredCurrentStatus > 3600*8
max_retries = 2
queue {{5}}
"""

    sh_file = """#!/bin/sh
unset PYTHONPATH
export X509_USER_PROXY=/afs/cern.ch/user/a/arnaik/x509up_u186038
"""

    redirectors = {
        "cms": "cmsxcache.crc.nd.edu",
        "fermi": "cmsxrootd.fnal.gov",
        "infn": "xrootd-cms.infn.it",
        "global": "cms-xrd-global.cern.ch",
    }


    # Make job submission directories!
    if not os.path.exists(os.path.join("/afs/cern.ch/user/a/arnaik/", ".event_mixing_submission")):
        os.mkdir(os.path.join("/afs/cern.ch/user/a/arnaik/", ".event_mixing_submission"))
    if not os.path.exists(os.path.join("/afs/cern.ch/user/a/arnaik/", ".event_mixing_submission", submission_dir)):
        os.mkdir(os.path.join("/afs/cern.ch/user/a/arnaik/", ".event_mixing_submission", submission_dir))
    if not os.path.exists(os.path.join("/afs/cern.ch/user/a/arnaik/", ".event_mixing_submission", submission_dir, "jobs")):
        os.mkdir(os.path.join("/afs/cern.ch/user/a/arnaik/", ".event_mixing_submission", submission_dir, "jobs"))
    if not os.path.exists(os.path.join("/afs/cern.ch/user/a/arnaik/", ".event_mixing_submission", submission_dir, "execs")):
        os.mkdir(os.path.join("/afs/cern.ch/user/a/arnaik/", ".event_mixing_submission", submission_dir, "execs"))
    if not os.path.exists(os.path.join("/afs/cern.ch/user/a/arnaik/", ".event_mixing_submission", submission_dir, "logs")):
        os.mkdir(os.path.join("/afs/cern.ch/user/a/arnaik/", ".event_mixing_submission", submission_dir, "logs"))

    # Confirm voms proxy before submitting jobs!
    #on_condorfe = subprocess.getstatusoutput("hostname")[1] == "condorfe.crc.nd.edu"
    on_condorfe = "cern.ch" in subprocess.getstatusoutput("hostname")[1]
    if not on_condorfe:
        stat, out = subprocess.getstatusoutput("voms-proxy-info -e -p")
        # stat is 0 the proxy is valid
        if stat != 0:
            raise RuntimeError("No valid proxy found. Please create one.")
        proxy = out

        # Set path of proxy as environment variables
        _x509_localpath = proxy
        _x509_path = os.environ["HOME"] + f'/.{_x509_localpath.split("/")[-1]}'
        os.system(f"cp {_x509_localpath} {_x509_path}")

    env_extra = [
        "export XRD_RUNFORKHANDLER=1",
        #f"export PYTHONPATH=$PYTHONPATH:/users/scastel2/miniforge3/envs/higgs-dna/bin/python",
        f"export PYTHONPATH=$HOME/.local/lib/python3.13/site-packages/:$PYTHONPATH",
    ]
    if not on_condorfe:
        env_extra.append(f"export X509_USER_PROXY={_x509_path}") 

    try:
        env_extra.append(f'export X509_CERT_DIR={os.environ["X509_CERT_DIR"]}')
    except KeyError:
        pass

    for cmd_str in env_extra:
        g_var = cmd_str.split()[-1]
        var, var_value = g_var.split("=")
        os.environ[var] = var_value


    # Split samples like HiggsDNA with sample JSONs
    sample_dict = None
    assert os.path.exists(os.path.abspath(os.path.join(cwd, args.file)))
    with open(os.path.abspath(os.path.join(cwd, args.file))) as f:
        sample_dict = json.load(f)

    # Use compiled version instead so there are no issues with loading multiple times!
    """
    # Write MakeClass macros here so we don't get an error when each job tries to overwrite Events.C/Event.h
    with open(os.path.join(cwd, "Events.C"), "w") as f:
        f.write(events_C)
    with open(os.path.join(cwd, "Events.h"), "w") as f:
        f.write(events_h)
    """

    # Write one sh and sub file for each era (key) in sample_dict. Each job will then get ONE queue number to run ONE of the if-statements in the sh file per file.
    for era in sample_dict:
        sub_path = f"condor_{era}.sub"
        sh_path = f"exec_{era}.sh"

        # One sh file for all jobs (files) in a given era
        with open(os.path.join("/afs/cern.ch/user/a/arnaik/", ".event_mixing_submission", submission_dir, "execs", sh_path), "w") as executable:
            executable.write(sh_file)
            executable.write("\n")

        # One sub file for all jobs (files) in a given era with queue of length N files * M cycles
        with open(os.path.join("/afs/cern.ch/user/a/arnaik/", ".event_mixing_submission", submission_dir, "jobs", sub_path), "w") as submission:
            submission.write(submit_file.format(
                sh_path,
                era,
                os.path.join("/afs/cern.ch/user/a/arnaik/", ".event_mixing_submission", submission_dir, "logs"),
                args.ram,
                args.queue,
                len(sample_dict[era])
            ))

        # Write each event mix file (job) to sh
        for idx, sample_file in enumerate(sample_dict[era]):
            # One job per file, cycles are handled in a loop in event mixing. Open files directly (no sample JSON needed)
            with open(os.path.join("/afs/cern.ch/user/a/arnaik/", ".event_mixing_submission", submission_dir, "execs", sh_path), "a") as executable:
                # Start if-statements
                executable.write(f"if [ $1 -eq {idx} ]; then\n")
                
                # Change redirector if specified
                sample_path = sample_file
                if args.xrootd is not None:
                    start_of_path = sample_file.find("///store")
                    sample_path = "root://" + redirectors[args.xrootd] + sample_file[start_of_path:]

                # Call python from /usr/bin/env to get my environment
                executable.write(f"/usr/bin/env python3 {os.path.abspath(os.path.join(cwd, 'do_event_mixing.py'))} {sample_path} {args.dump}/{era}/{sample_path.split('/')[-1]} {('--goldenJSON ' + os.path.abspath(args.goldenJSON)) if args.goldenJSON is not None else ''} --Ncycles {args.Ncycles} --cycles_offset {args.cycles_offset} {'--no-event-mixing' if args.no_event_mixing else ''} {'--debug' if args.debug else ''} {'--show-all' if args.show_all else ''} {'--noSkim' if args.noSkim else ''}".strip())
                executable.write("\n")
                executable.write(f"fi\n\n")

        if "root://" not in args.dump:
            if not os.path.exists(args.dump):
                Path(args.dump).mkdir(parents=True, exist_ok=True)

        # Set permissions of executable so condor can use it
        assert os.path.exists(os.path.join("/afs/cern.ch/user/a/arnaik/", ".event_mixing_submission", submission_dir, "execs", sh_path))
        os.system(f"chmod 775 {os.path.join('/afs/cern.ch/user/a/arnaik/', '.event_mixing_submission', submission_dir, 'execs', sh_path)}")

        # Submit jobs
        assert os.path.exists(os.path.join("/afs/cern.ch/user/a/arnaik/", ".event_mixing_submission", submission_dir, "jobs", sub_path))
        subprocess.run(["condor_submit", os.path.join("/afs/cern.ch/user/a/arnaik/", ".event_mixing_submission", submission_dir, "jobs", sub_path)])
