#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# fetching info from parquet files to mix photons
# --------------------------------------------------------------------------


def get_args():

    parser = argparse.ArgumentParser(
        description="Event mixing on post-baseprocessor parquet files."
    )

    parser.add_argument(
        "input",
        help="Merged parquet file from the base processor."
    )

    parser.add_argument(
        "outdir",
        help="Directory where mixed parquet will be saved."
    )

    parser.add_argument(
        "--nc",
        type=int,
        default=1,
        help="Number of mixing cycles."
    )

    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Starting cycle."
    )

    parser.add_argument(
        "--debug",
        action="store_true",
        help="Run only first 20 events."
    )

    return parser.parse_args()


def setup_logger(debug=False):

    level = logging.DEBUG if debug else logging.INFO

    logging.basicConfig(
        format="%(asctime)s : %(levelname)s : %(message)s",
        level=level,
    )

    return logging.getLogger(__name__)


def main():

    args = get_args()
    logger = setup_logger(args.debug)

    logger.info("Reading parquet file...")
    df = pd.read_parquet(args.input)
    logger.info(f"Loaded {len(df)} surviving events.")

    if args.debug:
        logger.warning("Running in DEBUG mode.")
        df = df.head(20).reset_index(drop=True)

    n_events = len(df)

    if n_events < 4:
        raise RuntimeError(
            "Need at least four surviving events for event mixing."
        )

    # ----------------------------------------------------------------
    # metadata bookkeeping
    # ----------------------------------------------------------------

    required = ["run", "luminosityBlock", "event"]

    for col in required:
        if col not in df.columns:
            raise RuntimeError(f"Missing column {col}")

    # fetching the pho{i} cols
    photon_columns = {}
    for pho in [1, 2, 3, 4]:
        photon_columns[pho] = [
            c for c in df.columns if c.startswith(f"pho{pho}_")
        ]

    # sanity check that pho2/3/4 columns match pho1's naming
    base_names = sorted(
        [c.replace("pho1_", "") for c in photon_columns[1]]
    )

    for pho in [2, 3, 4]:
        names = sorted(
            [c.replace(f"pho{pho}_", "") for c in photon_columns[pho]]
        )
        if names != base_names:
            raise RuntimeError(
                f"Photon {pho} columns do not match photon1."
            )

    # ----------------------------------------------------------------
    # Mixing algo
    # ----------------------------------------------------------------

    logger.info(f"Running {args.nc} mixing cycle(s)...")

    all_outputs = []
    skipped_cycles = 0
    max_cycle = n_events - 4

    indices = np.arange(n_events)

    for cycle in range(args.offset, args.offset + args.nc):

        logger.info(f"Processing cycle {cycle}")

        # Skip cycles that would produce duplicate photon combinations
        if cycle > max_cycle:
            logger.warning(
                f"Skipping cycle {cycle} "
                f"(would create duplicate photon combinations)"
            )
            skipped_cycles += 1
            continue

        logger.debug(f"Cycle {cycle}: creating mixed dataframe")

        mixed = df.copy(deep=True)
        mixed["mix_cycle"] = cycle

        # donor events (vectorized, no python-level row loop)
        j2 = (indices + 1 + cycle) % n_events
        j3 = (indices + 2 + cycle) % n_events
        j4 = (indices + 3 + cycle) % n_events

        for col in photon_columns[2]:
            mixed[col] = df.iloc[j2][col].to_numpy()
        for col in photon_columns[3]:
            mixed[col] = df.iloc[j3][col].to_numpy()
        for col in photon_columns[4]:
            mixed[col] = df.iloc[j4][col].to_numpy()

        all_outputs.append(mixed)

    # ----------------------------------------------------------------
    # combining all the event mixing cycles
    # ----------------------------------------------------------------

    if len(all_outputs) == 0:
        raise RuntimeError("No valid mixing cycles were produced.")

    logger.info("Combining all cycles...")

    mixed_df = pd.concat(all_outputs, ignore_index=True)
    mixed_df = mixed_df.reset_index(drop=True)

    logger.info(f"Final dataframe contains {len(mixed_df)} events.")

    # ----------------------------------------------------------------
    # writing output files
    # ----------------------------------------------------------------

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    input_name = Path(args.input).stem
    written_cycles = args.nc - skipped_cycles

    output_name = (
        f"{input_name}"
        f"_nc{written_cycles}"
        f"_offset{args.offset}.parquet"
    )

    output_path = outdir / output_name

    logger.info(f"Writing output to {output_path}")

    mixed_df.to_parquet(output_path, index=False)

    # ----------------------------------------------------------------
    # summary
    # ----------------------------------------------------------------

    logger.info("--------------------------------------------------")
    logger.info(f"Input events          : {len(df)}")
    logger.info(f"Output events         : {len(mixed_df)}")
    logger.info(f"Cycles requested      : {args.nc}")
    logger.info(f"Cycles written        : {written_cycles}")
    logger.info(f"Cycles skipped        : {skipped_cycles}")
    logger.info(f"Offset                : {args.offset}")
    logger.info(f"Output file           : {output_path}")
    logger.info("Finished successfully.")

    if skipped_cycles > 0:
        logger.warning(
            f"Skipped {skipped_cycles} cycle(s) "
            "to avoid duplicate events."
        )


if __name__ == "__main__":
    main()
