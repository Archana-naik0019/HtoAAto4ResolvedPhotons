import subprocess
import shutil
import os
import shlex
import argparse


if __name__ == "__main__":
    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser", conflict_handler="resolve")

    # Input arguments
    parser.add_argument("--cleanOnly", help="Only clean files without moving outputs to Combine folders. Defaults to false.", action="store_true", required=False)
    args = parser.parse_args()

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    for mass in range(15,65,5):
        print(f"~~~~   {mass} GeV   ~~~~")

        files = [
            f"higgsCombine{mass}GeV.HybridNew.mH125.root",
            f"higgsCombine{mass}GeV.HybridNew.mH125.quant0.500.root",
            f"higgsCombine{mass}GeV.HybridNew.mH125.quant0.025.root",
            f"higgsCombine{mass}GeV.HybridNew.mH125.quant0.160.root",
            f"higgsCombine{mass}GeV.HybridNew.mH125.quant0.840.root",
            f"higgsCombine{mass}GeV.HybridNew.mH125.quant0.975.root"
        ]

        # Compute limits
        print("Merging outputs...")

        if os.path.exists(f"higgsCombine{mass}GeV.HybridNew.mH125.merged.root"):
            os.remove(f"higgsCombine{mass}GeV.HybridNew.mH125.merged.root")
        for f in files:
            if not os.path.exists(f):
                files.remove(f)

        if len(files) != 0:
            # Merge outputs
            if not args.cleanOnly:
                command = f"hadd higgsCombine{mass}GeV.HybridNew.mH125.merged.root " + " ".join(files)
                subprocess.run(
                    shlex.split(command),
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    universal_newlines = True
                )

            # Move outputs to proper directory
            merged_path = f"{cwd}/../combine/{mass}_GeV/higgsCombine{mass}GeV.HybridNew.mH125.merged.root"
            if os.path.exists(merged_path):
                os.remove(merged_path)

            if not args.cleanOnly:
                shutil.move(f"higgsCombine{mass}GeV.HybridNew.mH125.merged.root", merged_path)

        # Removing logs and other outputs
        print("Removing output files...")
        for f in files:
            if os.path.exists(os.path.join(cwd,f)):
                os.remove(os.path.join(cwd, f))

        print("Removing extra condor outputs...")
        rm = os.listdir(cwd)
        for item in rm:
            #if (("obs" in item or "exp" in item or "68" in item or "95" in item) and \
            if "GeV" in item and \
            (".sub" in item or ".sh" in item or ".log" in item or ".err" in item or ".out" in item) or \
            "combine_logger.out" in item:
                os.remove(os.path.join(cwd, item))

        print()
