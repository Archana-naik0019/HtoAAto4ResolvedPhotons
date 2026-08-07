#!/usr/bin/env python3

"""
Loop over all ROOT files inside

input/
    2024C/
    2024D/
    ...

and invoke run_processor.py once per file.

Output structure:

output/
    2024C/
        file1_photons.parquet
        file1_cutflow.csv
        file2_photons.parquet
        file2_cutflow.csv
"""

import argparse
import glob
import os
import shutil
import subprocess


def get_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input-dir",
        required=True,
        help="Top-level directory containing era folders (2024C, 2024D, ...)",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        help="Top-level output directory",
    )

    parser.add_argument(
        "--metaconditions",
        required=True,
    )

    parser.add_argument(
        "--year",
        default="2024",
    )

    parser.add_argument(
        "--systematics",
        default="",
    )

    parser.add_argument(
        "--corrections",
        default="",
    )

    parser.add_argument(
        "--treename",
        default="Events",
    )

    parser.add_argument(
        "--doFlow_corrections",
        action="store_true",
    )

    parser.add_argument(
        "--use_photonid_mva",
        action="store_true",
    )

    parser.add_argument(
        "--debug",
        action="store_true",
    )

    return parser


def main():

    args = get_parser().parse_args()

    eras = sorted(
        d for d in os.listdir(args.input_dir)
        if os.path.isdir(os.path.join(args.input_dir, d))
    )

    for era in eras:

        input_era = os.path.join(args.input_dir, era)
        output_era = os.path.join(args.output_dir, era)

        os.makedirs(output_era, exist_ok=True)

        root_files = sorted(glob.glob(os.path.join(input_era, "*.root")))

        print(f"\nProcessing {era} ({len(root_files)} files)")

        for root_file in root_files:

            basename = os.path.splitext(os.path.basename(root_file))[0]

            tmp_output = os.path.join(output_era, "__tmp__")

            if os.path.exists(tmp_output):
                shutil.rmtree(tmp_output)

            os.makedirs(tmp_output)

            cmd = [
                "python",
                "run_processor.py",
                "--input", root_file,
                "--output", tmp_output,
                "--dataset", era,
                "--year", args.year,
                "--metaconditions", args.metaconditions,
                "--systematics", args.systematics,
                "--corrections", args.corrections,
                "--treename", args.treename,
            ]

            if args.doFlow_corrections:
                cmd.append("--doFlow_corrections")

            if args.use_photonid_mva:
                cmd.append("--use_photonid_mva")

            if args.debug:
                cmd.append("--debug")

            print(f"  -> {os.path.basename(root_file)}")

            subprocess.run(cmd, check=True)

            shutil.move(
                os.path.join(tmp_output, "photons.parquet"),
                os.path.join(output_era, f"{basename}_photons.parquet"),
            )

            shutil.move(
                os.path.join(tmp_output, "cutflow.csv"),
                os.path.join(output_era, f"{basename}_cutflow.csv"),
            )

            shutil.rmtree(tmp_output)


if __name__ == "__main__":
    main()
