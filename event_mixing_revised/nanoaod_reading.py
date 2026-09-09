import re
import awkward as ak
import numpy as np
import uproot


def group_branches(branches: list) -> tuple:
    """
    Split a flat list of NanoAOD branch names into:
      - collections
      - flat (no counter associated)
    """

    def is_counter(b):
        return b.startswith("n") and len(b) > 1 and b[1].isupper()

    counters = sorted(b for b in branches if is_counter(b))
    collections = {}
    claimed = set()

    for counter in counters:
        coll = counter[1:]
        prefix = f"{coll}_"
        members = [b for b in branches if b.startswith(prefix)]
        if members:
            collections[coll] = members
            claimed.add(counter)
            claimed.update(members)

    flat = [b for b in branches if b not in claimed]
    return collections, flat


def read_nanoaod_events(
    filepaths,
    treename: str = "Events",
    entry_stop: int = None,
    branch_filter=None,
) -> tuple:
    """
    Read one or more NanoAOD ROOT files into a single nested awkward Array,
    with jagged collections (Photons) grouped

    Parameters------------------------------------------------
    filepaths : str or list[str]
    treename : str, default "Events"
    ----------------------------------------------------------

    Returns---------------------------------------------------
    events : ak.Array
    collections : dict[str, list[str]]
        As returned by group_branches, useful for logging/debugging.
    flat : list[str]
    dtypes : dict[str, np.dtype or None]
        Per-branch numpy dtype as read from the input file(s), captured
        *before* any zipping/reshuffling. Pass this to restore_dtypes()
        right before writing.
    """
    if isinstance(filepaths, str):
        filepaths = [filepaths]

    with uproot.open(f"{filepaths[0]}:{treename}") as t0:
        all_branches = t0.keys()

    if branch_filter is not None:
        all_branches = [b for b in all_branches if branch_filter(b)]

    collections, flat = group_branches(all_branches)

    files = {fp: treename for fp in filepaths}
    arrays = uproot.concatenate(
        files, filter_name=all_branches, entry_stop=entry_stop, library="ak"
    )

    dtypes = {b: _infer_dtype(arrays[b]) for b in all_branches}

    fields = {}
    for b in flat:
        fields[b] = arrays[b]

    for coll, members in collections.items():
        sub = {}
        for m in members:
            field_name = m[len(coll) + 1:]
            sub[field_name] = arrays[m]
        fields[coll] = ak.zip(sub, depth_limit=1)

    events = ak.zip(fields, depth_limit=1)
    return events, collections, flat, dtypes


def _infer_dtype(array: ak.Array, sample_size: int = 2000):
    """
    Finds dtypes
    """
    sample = array[: min(sample_size, len(array))]
    flat = ak.flatten(sample, axis=None)
    if len(flat) == 0:
        flat = ak.flatten(array, axis=None)
    if len(flat) == 0:
        return None
    try:
        return ak.to_numpy(flat).dtype
    except Exception:
        return None


def restore_dtypes(
    events: ak.Array, dtypes: dict, collections: dict, flat: list
) -> ak.Array:
    """
    Cast every branch that existed in the *input* schema back to the numpy
    dtype it had on read, undoing any int64/float64 upcasting introduced
    by intermediate awkward operations (concatenate, indexing etc. can silently promote dtypes)

    called right before write_root_events() so the output file matches
    the input file's schema.
    """
    out = events
    existing = set(ak.fields(out))

    for b in flat:
        dt = dtypes.get(b)
        if dt is not None and b in existing:
            out = ak.with_field(out, ak.values_astype(out[b], dt), b)

    for coll, members in collections.items():
        if coll not in existing:
            continue
        sub_fields = set(ak.fields(out[coll]))
        new_sub = {}
        changed = False
        for m in members:
            field_name = m[len(coll) + 1:]
            if field_name not in sub_fields:
                continue
            dt = dtypes.get(m)
            if dt is not None:
                new_sub[field_name] = ak.values_astype(out[coll][field_name], dt)
                changed = True
            else:
                new_sub[field_name] = out[coll][field_name]
        if changed:
            out = ak.with_field(out, ak.zip(new_sub, depth_limit=1), coll)

    return out


def drop_field(events: ak.Array, name: str) -> ak.Array:
    """Remove a top-level field (the added 'mix_cycle') if present."""
    if name not in ak.fields(events):
        return events
    remaining = {f: events[f] for f in ak.fields(events) if f != name}
    return ak.zip(remaining, depth_limit=1)


def skim_events(
    events: ak.Array,
    n_photons: int = 4,
    hlt_branch: str = "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId",
    no_skim: bool = True,
) -> ak.Array:
    """
    Equivalent of TreeHelper.SkimTree: require >= n_photons photons, and
    (unless no_skim) that the given HLT bit fired.
    """
    if no_skim:
        return events

    has_min_photons = ak.num(events.Photon.pt, axis=1) >= n_photons

    if hlt_branch is None:
        mask = has_min_photons
    else:
        if hlt_branch not in ak.fields(events):
            raise KeyError(
                f"HLT branch '{hlt_branch}' not found in events. "
                "Pass --noSkim or a different --hlt-branch."
            )
        mask = has_min_photons & (events[hlt_branch] == 1)

    return events[mask]


def truncate_photons(events: ak.Array, n_photons: int = 4) -> ak.Array:
    
    photon_fields = ak.fields(events.Photon)
    trunc = {f: events.Photon[f][:, :n_photons] for f in photon_fields}
    trunc_zip = ak.zip(trunc, depth_limit=1)
    trunc_zip = ak.from_regular(trunc_zip, axis=1)  # keeping it jagged

    return ak.with_field(events, trunc_zip, "Photon")

#--------------this version of write_root_events crashes at systematics in HiggsDNA-------------
'''
def write_root_events(events: ak.Array, outpath: str, treename: str = "Events") -> None:
    """
    Writes an Awkward events array out explicitly as a legacy ROOT TTree.
    """
    out_dict = {}

    for field in ak.fields(events):
        col = events[field]
        subfields = ak.fields(col)

        # 1. Jagged 2D Collections (Electron, Photon, Jet)
        if len(subfields) > 0:
            first_subfield = subfields[0]
            counter_name = f"n{field}"

            if counter_name not in ak.fields(events) and counter_name not in out_dict:
                out_dict[counter_name] = ak.to_numpy(ak.num(col[first_subfield], axis=1)).astype(np.int32)

            for subfield in subfields:
                out_dict[f"{field}_{subfield}"] = ak.to_packed(col[subfield])

        # 2. 1D Scalar Branches
        else:
            out_dict[field] = ak.to_numpy(col)

    # Force classic TTree serialization
    with uproot.recreate(outpath) as f:
        # Build branch type mapping to force TTree backend
        branch_types = {}
        for k, v in out_dict.items():
            if isinstance(v, np.ndarray):
                branch_types[k] = v.dtype
            else:
                branch_types[k] = v.type

        # mktree creates the TTree metadata key explicitly
        f.mktree(treename, branch_types)
        f[treename].extend(out_dict)
'''
#---------------------------------------------------------------------------------------
def write_root_events(events: ak.Array, outpath: str, treename: str = "Events") -> None:
    out_dict = {}

    for field in ak.fields(events):
        col = events[field]
        subfields = ak.fields(col)

        if subfields:
            # Jagged collection (Photon, Electron, Jet, ...).
            # Computes ONE shared per-event count and derives every sub-branch from the SAME col, so uproot recognizes them as sharing one counter instead of minting one per field.
            counts = ak.num(col[subfields[0]], axis=1)
            out_dict[f"n{field}"] = ak.to_numpy(counts).astype(np.int32)
            for subfield in subfields:
                out_dict[f"{field}_{subfield}"] = col[subfield]
        else:
            out_dict[field] = ak.to_numpy(col) if not ak.fields(col) else col

    with uproot.recreate(outpath) as f:
        # Force classic TTree serialization (not RNTuple). Without forcing this the output assumed RNTuple form.
        # systematics-variation arrays require this explicitly (rest of the framework works fine even without this)*****************************
        branch_types = {}
        for k, v in out_dict.items():
            if isinstance(v, np.ndarray):
                branch_types[k] = v.dtype
            else:
                branch_types[k] = ak.Array(v).type
        f.mktree(treename, branch_types)
        #****************************************************************************************************************************************
        
        f[treename].extend(out_dict)
#-------------------------------------------------------------------------------------------


def write_parquet_events(events: ak.Array, outpath: str) -> None:
    ak.to_parquet(events, outpath)


def make_output_path(outpath: str, n_cycles: int, cycles_offset: int, ext: str) -> str:
    base = re.sub(r"\.(root|parquet)$", "", outpath)
    return f"{base}_nc{n_cycles}_offset{cycles_offset}.{ext}"
