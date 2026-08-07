#!/usr/bin/env python3

import os
import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

#--------------------

#==============================fetching info from parquet files to mix photons=======================

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

#---------------------------

def setup_logger(debug=False):

    level = logging.DEBUG if debug else logging.INFO

    logging.basicConfig(
        format="%(asctime)s : %(levelname)s : %(message)s",
        level=level,
    )

    return logging.getLogger(__name__)

#----------------------------------

def main():

    args = get_args()

    logger = setup_logger(args.debug)

    logger.info("Reading parquet file...")

    df = pd.read_parquet(args.input)

    logger.info(f"Loaded {len(df)} surviving events.")

#--------------------------------------

# this section is redundant..may be I will remove it after checks
    if args.debug:

        logger.warning("Running in DEBUG mode.")

        df = df.head(20)

#------------------------------------------

    N = len(df)

    if N < 4:
        raise RuntimeError(
            "Need at least four surviving events for event mixing."
        )
#-----------------------------------------------

# metadata bookkeeping

required = [
    "run",
    "luminosityBlock",
    "event"
]

for col in required:

    if col not in df.columns:

        raise RuntimeError(f"Missing column {col}")

#------------------------------------------------

# fetching the pho{i} cols

photon_columns = {}

for pho in [1,2,3,4]:

    photon_columns[pho] = [
        c for c in df.columns
        if c.startswith(f"pho{pho}_")
    ]

#---------------------------------------------------

# a sanity check, might me removed later

base_names = sorted(
    [c.replace("pho1_", "") for c in photon_columns[1]]
)

for pho in [2,3,4]:

    names = sorted(
        [c.replace(f"pho{pho}_", "") for c in photon_columns[pho]]
    )

    if names != base_names:

        raise RuntimeError(
            f"Photon {pho} columns do not match photon1."
        )

#-------------------------------------------------

#===========================================Mixing algo===============================================

mixed_dfs = []

logger.info(f"Running {args.nc} mixing cycle(s)...")

N = len(df)

# looping over cycles
for cycle in range(args.offset, args.offset + args.nc):

    logger.info(f"Processing cycle {cycle}")
    logger.debug(f"Cycle {cycle}: creating mixed dataframe")

    # copying the dataframe for this cycle
    mixed = df.copy(deep=True)
    mixed["mix_cycle"] = cycle
    indices = np.arange(N)

    # donor events defined
    j2 = (indices + 1 + cycle) % N
    j3 = (indices + 2 + cycle) % N
    j4 = (indices + 3 + cycle) % N

#------------------------------------------------------------------

for col in photon_columns[2]:

    mixed[col] = df.iloc[j2][col].to_numpy()
for col in photon_columns[3]:

    mixed[col] = df.iloc[j3][col].to_numpy()
for col in photon_columns[4]:

    mixed[col] = df.iloc[j4][col].to_numpy()

#-----------------------------------------------------------------

# needs rechecking

skipped_cycles = 0

all_outputs = []

for cycle in range(args.cycles_offset,
                   args.cycles_offset + args.Ncycles):

    logger.info(f"Processing cycle {cycle}")

    # Skip cycles that would produce duplicate events
    max_cycle = len(df) - 4

    if cycle > max_cycle:

        logger.warning(
            f"Skipping cycle {cycle} "
            f"(would create duplicate photon combinations)"
        )

        skipped_cycles += 1
        continue

    # Start with an exact copy of the dataframe
    mixed = df.copy()

    logger.debug(f"Cycle {cycle}: creating mixed dataframe")

    # Keep track of which mixing cycle produced this row
    mixed["mix_cycle"] = cycle

    n_events = len(df)

    # ------------------------------------------------------------------
    # Replace photons 2,3,4
    # ------------------------------------------------------------------

    for i in range(n_events):

        idx2 = (i + cycle + 1) % n_events
        idx3 = (i + cycle + 2) % n_events
        idx4 = (i + cycle + 3) % n_events

        # Photon 2
        for col in photon_columns[2]:
            mixed.at[i, col] = df.at[idx2, col]

        # Photon 3
        for col in photon_columns[3]:
            mixed.at[i, col] = df.at[idx3, col]

        # Photon 4
        for col in photon_columns[4]:
            mixed.at[i, col] = df.at[idx4, col]

    all_outputs.append(mixed)

#-----------------------------------------------------------------------

# combining all the event mixing cycles

if len(all_outputs) == 0:
    raise RuntimeError("No valid mixing cycles were produced.")

logger.info("Combining all cycles...")

mixed_df = pd.concat(
    all_outputs,
    ignore_index=True
)

logger.info(
    f"Final dataframe contains {len(mixed_df)} events."
)

#-----------------------------------------------------------------------

mixed_df = mixed_df.reset_index(drop=True)

#----------------------------------------------------------------

#writing output files

outdir = Path(args.outpath)
outdir.mkdir(parents=True, exist_ok=True)

# Input filename without extension
input_name = Path(args.file).stem
written_cycles = args.Ncycles - skipped_cycles

output_name = (
    f"{input_name}"
    f"_nc{written_cycles}"
    f"_offset{args.cycles_offset}.parquet"
)

output_path = outdir / output_name

logger.info(f"Writing output to {output_path}")

mixed_df.to_parquet(
    output_path,
    index=False
)

#------------------------------------------------------------------
# summary

logger.info("--------------------------------------------------")
logger.info(f"Input events          : {len(df)}")
logger.info(f"Output events         : {len(mixed_df)}")

written_cycles = args.Ncycles - skipped_cycles

logger.info(f"Cycles requested : {args.Ncycles}")
logger.info(f"Cycles written   : {written_cycles}")
logger.info(f"Cycles skipped   : {skipped_cycles}")
logger.info(f"Cycles written        : {args.Ncycles}")
logger.info(f"Offset                : {args.cycles_offset}")
logger.info(f"Output file           : {output_path}")
logger.info("Finished successfully.")

# warning if any cycles had to be skipped
if skipped_cycles > 0:
    logger.warning(
        f"Skipped {skipped_cycles} cycle(s) "
        "to avoid duplicate events."
    )
