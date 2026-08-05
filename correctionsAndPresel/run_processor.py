#!/usr/bin/env python
"""
driver for HggBaseProcessor (minimal_base_2.py).
Loads a single ROOT file, builds the events array, and runs
process() directly - no dask/parsl/coffea.Runner needed since
this processor writes its own parquet/csv output rather than
returning accumulatable histograms.

commandline to be used: python3 run_processor.py   --input /eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/event_mix_2024Full_test1_nc2_offset0/2024C/dbc7f6a4-74fc-40e2-81ff-0e70d2eabe63_nc2_offset0.root   --output /eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/correctionsAndPresel/corr_applied_testing1   --year 2024   --metaconditions metaconditions.json   --corrections Scale_IJazZ
"""

import argparse
import json
import logging
import os

from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

from minimal_base_2 import HggBaseProcessor


def get_parser():
    parser = argparse.ArgumentParser(description="Run HggBaseProcessor on a single ROOT file")
    parser.add_argument(
        "--input", required=True, type=str,
        help="Path to input ROOT file (local or root:// / xrootd path)",
    )
    parser.add_argument(
        "--output", required=True, type=str,
        help="Output directory for photons.parquet and cutflow.csv",
    )
    parser.add_argument(
        "--dataset", type=str, default=None,
        help="Dataset name (used as key into --year/--systematics/--corrections). "
             "Defaults to the input file's parent directory name.",
    )
    parser.add_argument(
        "--year", type=str, default="2024",
        help="Year/era string, e.g. 2024, 2022EE, 2022postEE",
    )
    parser.add_argument(
        "--metaconditions", required=True, type=str,
        help="Path to metaconditions JSON file",
    )
    parser.add_argument(
        "--systematics", type=str, default="",
        help="Comma-separated list of systematic names to apply (default: none)",
    )
    parser.add_argument(
        "--corrections", type=str, default="",
        help="Comma-separated list of correction names to apply (default: none)",
    )
    parser.add_argument(
        "--doFlow_corrections", action="store_true",
        help="Enable normalizing-flow corrections (MC only)",
    )
    parser.add_argument(
        "--use_photonid_mva", action="store_true",
        help="Recompute photon ID MVA on the fly (skip for data-only runs that don't need it)",
    )
    parser.add_argument(
        "--treename", type=str, default="Events",
        help="TTree name inside the ROOT file (default: Events)",
    )
    parser.add_argument(
        "--debug", action="store_true",
        help="Enable debug-level logging",
    )
    return parser


def main():
    args = get_parser().parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    logger = logging.getLogger(__name__)

    # Derive dataset name if not given explicitly
    if args.dataset is None:
        # e.g. .../event_mix_2024Full_test1_nc2_offset0/2024C/file.root -> "2024C"
        dataset_name = os.path.basename(os.path.dirname(args.input))
    else:
        dataset_name = args.dataset
    logger.info(f"Using dataset name: {dataset_name}")

    with open(args.metaconditions) as f:
        metaconditions = json.load(f)

    systematics_list = [s for s in args.systematics.split(",") if s]
    corrections_list = [c for c in args.corrections.split(",") if c]

    systematics = {dataset_name: systematics_list} if systematics_list else {}
    corrections = {dataset_name: corrections_list} if corrections_list else {}
    year = {dataset_name: [args.year]}

    logger.info(f"Systematics: {systematics}")
    logger.info(f"Corrections: {corrections}")
    logger.info(f"Year: {year}")
    logger.info(f"use_photonid_mva: {args.use_photonid_mva}")
    logger.info(f"doFlow_corrections: {args.doFlow_corrections}")

    os.makedirs(args.output, exist_ok=True)

    processor_instance = HggBaseProcessor(
        metaconditions=metaconditions,
        systematics=systematics,
        corrections=corrections,
        year=year,
        output_location=args.output,
        doFlow_corrections=args.doFlow_corrections,
        use_photonid_mva=args.use_photonid_mva,
    )

    logger.info(f"Loading events from: {args.input}")
    events = NanoEventsFactory.from_root(
        {args.input: args.treename},
        schemaclass=NanoAODSchema,
        metadata={"dataset": dataset_name},
    ).events()

    logger.info("Running process()...")
    result = processor_instance.process(events)

    logger.info(f"Done. Result summary: {result}")


if __name__ == "__main__":
    main()
