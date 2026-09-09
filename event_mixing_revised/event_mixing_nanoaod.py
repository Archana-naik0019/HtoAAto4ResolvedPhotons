import awkward as ak
import numpy as np


def mix_photons_nanoaod(
    events: ak.Array,
    n_cycles: int = 20,
    cycle_offset: int = 0,
    n_photons: int = 4,
) -> ak.Array:
    """
    Performs event mixing on preselected events using Awkward Arrays.
    Preserves all NanoAOD fields attached to the Photon record collection,
    as well as every other branch/collection on `events` (jets, electrons,
    HLT bits, etc.), which are copied unchanged per cycle.
    """
    n_events = len(events)
    if n_events < n_photons:
        raise ValueError(
            f"Need at least {n_photons} surviving events to perform event mixing "
            f"(got {n_events})."
        )

    base_indices = np.arange(n_events)
    photon_fields = ak.fields(events.Photon)

    mixed_cycle_list = []
    skipped_cycles = 0

    for c in range(cycle_offset, cycle_offset + n_cycles):
        # Avoid a mixed photon landing back on its own event.
        if c > (n_events - n_photons):
            skipped_cycles += 1
            continue

        # index arrays for photon slots 1..n_photons-1; slot 0 is unshifted
        shot_indices = [base_indices]
        for k in range(1, n_photons):
            shot_indices.append((base_indices + k + c) % n_events)

        mixed_photon_dict = {}
        for field in photon_fields:
            vals = [events.Photon[field][shot_indices[k], k] for k in range(n_photons)]
            mixed_photon_dict[field] = ak.concatenate(
                [v[:, None] for v in vals], axis=1
            )

        mixed_photons = ak.zip(mixed_photon_dict, depth_limit=1)
        mixed_photons = ak.from_regular(mixed_photons, axis=1) # added to avoid a crash
        cycle_events = ak.with_field(events, mixed_photons, "Photon")
        cycle_events = ak.with_field(
            cycle_events, ak.full_like(events.event, c, dtype=np.int32), "mix_cycle"
        )
        mixed_cycle_list.append(cycle_events)

    if not mixed_cycle_list:
        raise ValueError(
            "No valid mixing cycles could be produced - n_cycles/cycle_offset "
            "too large relative to the number of surviving events."
        )

    return ak.concatenate(mixed_cycle_list, axis=0), skipped_cycles
