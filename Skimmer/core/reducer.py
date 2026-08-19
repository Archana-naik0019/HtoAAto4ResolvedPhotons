import awkward as ak
from core.config import COLLECTIONS, SCALARS, SCALARS_MC, WEIGHTS
from core.event_store import EventStore
from core.reader import NanoReader
from selection.electron import electron_mask
from selection.event import event_mask
from selection.jet import jet_mask
from selection.muon import muon_mask
from selection.photon import photon_mask


class NanoReducer:

    def __init__(
        self,
        input_file,
        jet_selection=True,
        electron_selection=True,
        muon_selection=True,
        photon_selection=True,
        event_selection=True,
    ):
        self.reader = NanoReader(input_file)
        self.jet_selection = jet_selection
        self.electron_selection = electron_selection
        self.muon_selection = muon_selection
        self.photon_selection = photon_selection
        self.event_selection = event_selection

    def run(self):

        collections = {}
        store = EventStore()

        #
        # Dynamically detect if this file is MC or Data
        #
        try:
            genWeight = self.reader.read_scalar("genWeight")
            is_mc = True
        except Exception:
            genWeight = None
            is_mc = False

        if is_mc:
            original_index = ak.local_index(genWeight)
            store.add_temp("genWeight_original", genWeight)
            store.add_temp("original_index", original_index)

            store.add_metadata("sum_genw_presel", float(ak.sum(genWeight)))
            store.add_metadata("n_events_presel", len(genWeight))
        else:
            # For Data: use 'event' or 'run' branch to extract the total event count/index
            run_num = self.reader.read_scalar("run")
            original_index = ak.local_index(run_num)

            store.add_temp("original_index", original_index)
            store.add_metadata("n_events_presel", len(run_num))

        #
        # Read collections and apply object selections
        #
        for collection in COLLECTIONS:
            print(f"Reading {collection}")
            obj = self.reader.read(collection)

            if collection == "Jet" and self.jet_selection:
                obj = obj[jet_mask(obj)]
            elif collection == "Electron" and self.electron_selection:
                obj = obj[electron_mask(obj)]
            elif collection == "Muon" and self.muon_selection:
                obj = obj[muon_mask(obj)]
            elif collection == "Photon" and self.photon_selection:
                obj = obj[photon_mask(obj)]

            collections[collection] = obj

        # Event selection
        if self.event_selection:
            mask = event_mask(collections)
        else:
            mask = ak.ones_like(original_index, dtype=bool)

        store.add_scalar("__original_index__", original_index[mask])

        # Save collections
        for name, obj in collections.items():
            store.add_collection(name, obj[mask])

        # Save standard scalars (Common to Data and MC)
        for branch in SCALARS:
            print(f"Reading {branch}")
            try:
                value = self.reader.read_scalar(branch)
                store.add_scalar(branch, value[mask])
            except (KeyError, Exception) as e:
                print(f"Skipping missing scalar '{branch}' for this dataset.")

        # Save MC-only scalars and weights if present
        if is_mc:
            for branch in SCALARS_MC:
                print(f"Reading {branch}")
                try:
                    value = self.reader.read_scalar(branch)
                    store.add_scalar(branch, value[mask])
                except Exception:
                    pass

            for branch in WEIGHTS:
                print(f"Reading {branch}")
                try:
                    value = self.reader.read_weight(branch)
                    store.add_weight(branch, value[mask])
                except Exception:
                    pass

        return store
