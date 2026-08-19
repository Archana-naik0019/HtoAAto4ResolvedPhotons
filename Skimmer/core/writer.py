import os

import awkward as ak
import numpy as np
import uproot

from core.config import MAX_EVENTS_PER_FILE


class NanoWriter:

    def __init__(self, filename):
        self.filename = filename

    def write(self, store):

        # Check if genWeight_original exists in temp store (MC vs Data)
        has_genweight = "genWeight_original" in store.temp
        if has_genweight:
            genWeight_original = store.temp["genWeight_original"]

        # Original indices of surviving events
        original_index = ak.to_numpy(store.scalars["__original_index__"])

        n_skim = len(original_index)

        #
        # Determine skimmed split points
        #
        if n_skim <= MAX_EVENTS_PER_FILE:
            split_points = []
        else:
            split_points = np.arange(
                MAX_EVENTS_PER_FILE,
                n_skim,
                MAX_EVENTS_PER_FILE,
            )

        starts = np.concatenate((np.array([0]), split_points))
        stops = np.concatenate((split_points, np.array([n_skim])))

        #
        # Original NanoAOD boundaries
        #
        original_start = 0

        for i, (start, stop) in enumerate(zip(starts, stops)):

            #
            # Last file
            #
            if stop == n_skim:
                # If we have genWeight, use its total length, otherwise infer from scalars/temp or max original index
                if has_genweight:
                    original_stop = len(genWeight_original)
                else:
                    # Fallback for data: if n_skim > 0, highest original index + 1
                    original_stop = (
                        int(original_index[-1]) + 1 if n_skim > 0 else 0
                    )
            else:
                # Original event corresponding to last skimmed event
                original_stop = original_index[stop - 1] + 1

            #
            # Slice branches
            #
            branches = {}

            for name, array in store.scalars.items():

                if name == "__original_index__":
                    continue

                branches[name] = array[start:stop]

            for name, array in store.weights.items():
                branches[name] = array[start:stop]

            for name, array in store.collections.items():
                branches[name] = array[start:stop]

            #
            # Metadata
            #
            metadata = {
                "n_events_presel": np.array(
                    [original_stop - original_start], dtype=np.int64
                ),
            }

            if has_genweight:
                metadata["sum_genw_presel"] = np.array(
                    [
                        float(
                            ak.sum(
                                genWeight_original[
                                    original_start:original_stop
                                ]
                            )
                        )
                    ],
                    dtype=np.float64,
                )

            #
            # Output filename
            #
            if len(starts) == 1:
                outfile = self.filename
            else:
                base, ext = os.path.splitext(self.filename)
                outfile = f"{base}_{i:03d}{ext}"

            #
            # Write ROOT file
            #
            with uproot.recreate(
                outfile,
                compression=uproot.ZSTD(9),
            ) as fout:

                events = fout.mktree(
                    "Events",
                    {
                        name: array.type
                        for name, array in branches.items()
                    },
                )

                events.extend(branches)

                meta = fout.mktree(
                    "Metadata",
                    {
                        key: value.dtype
                        for key, value in metadata.items()
                    },
                )

                meta.extend(metadata)

            print(
                f"Wrote {outfile}"
                f" | skimmed events = {stop-start}"
                f" | original events = {original_stop-original_start}"
            )

            #
            # Next chunk starts here
            #
            original_start = original_stop
