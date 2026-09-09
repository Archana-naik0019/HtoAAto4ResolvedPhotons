#!/usr/bin/env python

import argparse
import logging
import os
import sys

import awkward as ak

# imports from nanoaod_reading.py
from nanoaod_reading import (
    read_nanoaod_events,
    skim_events,
    truncate_photons,
    restore_dtypes,
    drop_field,
    write_root_events,
    write_parquet_events,
    make_output_path,
)
from event_mixing_nanoaod import mix_photons_nanoaod


class SplitArgs(argparse.Action):
    """Splits comma-separated strings into lists for argparse."""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(","))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=str, action=SplitArgs,
                         help="File(s) to process. Comma-separated for multiple paths.")
    parser.add_argument("outpath", type=str, help="Output path (extension is replaced).")
    parser.add_argument("-nc", "--Ncycles", type=int, default=1, required=True,
                         help="Number of event-mixing cycles to run.")
    parser.add_argument("-offset", "--cycles_offset", type=int, default=0, required=False,
                         help="Starting offset for mixing cycles. Defaults to 0.")
    parser.add_argument("--n-photons", type=int, default=4,
                         help="Number of leading photons to keep/mix per event. Defaults to 4.")
    parser.add_argument("--hlt-branch", type=str,
                         default="HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId",
                         help="HLT branch required to pass the skim (ignored with --noSkim).")
    parser.add_argument("--no-event-mixing", action="store_true",
                         help="Skip mixing; just skim and truncate to n_photons.")
    parser.add_argument("--noSkim", action="store_true",
                         help="Skip the HLT requirement; only apply the photon-count cut.")
    parser.add_argument("--strict-schema", action="store_true",
                         help="Drop added fields (like mix_cycle) so schema matches input 100%.")
    parser.add_argument("--outtype", choices=["root", "parquet", "both"], default="root",
                         help="Output format(s) to write. Defaults to root.")
    parser.add_argument("--treename", type=str, default="Events",
                         help="Input/output TTree name. Defaults to 'Events'.")
    parser.add_argument("--debug", action="store_true",
                         help="Run in debug mode. Limits to 20 skimmed events!")
    args = parser.parse_args()

    assert args.Ncycles >= 1

    logger = logging.getLogger()
    logging.basicConfig(
        stream=sys.stdout,
        format="%(asctime)s: %(levelname)s: %(message)s",
        level=logging.INFO if not args.debug else logging.DEBUG,
    )

    for f in args.file:
        if "root://" not in f:
            assert os.path.exists(f), f"Input file not found: {f}"

    logger.info(f"Processing file(s): {', '.join(os.path.basename(f) for f in args.file)}")

    logger.info("Reading input file(s)...")
    # Capturing all 4 return values (including original branch dtypes)
    events, collections, flat, dtypes = read_nanoaod_events(args.file, treename=args.treename)
    logger.info(f"Loaded {len(events)} raw events with {len(collections)} collections "
                f"and {len(flat)} flat branches.")

    logger.info("Applying skim (photon count%s)..." %
                (" + HLT" if not args.noSkim else ""))
    events = skim_events(
        events,
        n_photons=args.n_photons,
        hlt_branch=args.hlt_branch,
        no_skim=args.noSkim,
    )
    n_skimmed = len(events)
    logger.info(f"Skimmed down to {n_skimmed} events.")

    if n_skimmed == 0:
        logger.info("0 events survived the skim. Exiting.")
        sys.exit(0)

    if args.debug:
        events = events[: min(20, n_skimmed)]
        logger.debug(f"Debug mode: limited to {len(events)} events.")

    if args.no_event_mixing:
        logger.warning("Specified NO event mixing. Truncating to n_photons only.")
        out_events = truncate_photons(events, n_photons=args.n_photons)
        n_cycles_used, skipped_cycles = 1, 0
    else:
        logger.info(f"Running event mixing: {args.Ncycles} cycle(s) from offset "
                    f"{args.cycles_offset}, keeping {args.n_photons} photons/event.")
        out_events, skipped_cycles = mix_photons_nanoaod(
            events,
            n_cycles=args.Ncycles,
            cycle_offset=args.cycles_offset,
            n_photons=args.n_photons,
        )
        n_cycles_used = args.Ncycles - skipped_cycles
        if skipped_cycles > 0:
            logger.warning(f"Required to skip {skipped_cycles} cycle(s) to avoid "
                            "duplicate/self-mixed photons.")

    logger.info(f"Output event count: {len(out_events)}.")

    if len(out_events) == 0:
        logger.info("Nothing to write (0 output events). Exiting.")
        sys.exit(0)

    # Clean up and restoring exact input data types prior to writing
    if args.strict_schema:
        out_events = drop_field(out_events, "mix_cycle")
        
    out_events = restore_dtypes(out_events, dtypes, collections, flat)

    if args.outtype in ("root", "both"):
        root_out = make_output_path(args.outpath, n_cycles_used, args.cycles_offset, "root")
        os.makedirs(os.path.dirname(os.path.abspath(root_out)) or ".", exist_ok=True)
        logger.info(f"Writing ROOT output to {root_out}")
        write_root_events(out_events, root_out, treename=args.treename)
        assert os.path.exists(root_out) and os.stat(root_out).st_size != 0, \
            "ROOT output writing failed!"

    if args.outtype in ("parquet", "both"):
        pq_out = make_output_path(args.outpath, n_cycles_used, args.cycles_offset, "parquet")
        os.makedirs(os.path.dirname(os.path.abspath(pq_out)) or ".", exist_ok=True)
        logger.info(f"Writing parquet output to {pq_out}")
        write_parquet_events(out_events, pq_out)
        assert os.path.exists(pq_out) and os.stat(pq_out).st_size != 0, \
            "Parquet output writing failed!"

    logger.info(f"Finished processing {len(args.file)} file(s): "
                f"{n_skimmed} skimmed -> {len(out_events)} output events "
                f"from {n_cycles_used} cycle(s).")


if __name__ == "__main__":
    main()
