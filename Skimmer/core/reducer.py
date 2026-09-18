import awkward as ak
from core.config import COLLECTIONS, SCALARS, SCALARS_MC, WEIGHTS, KEEP_FIELDS
from core.event_store import EventStore
from core.reader import NanoReader
from selection.electron import electron_mask
from selection.event import event_mask
#from selection.photon import photon_mask
from selection.photon import h4g_preselection_mask
from selection.photon import (
    GAP_BARREL_ETA,
    GAP_ENDCAP_ETA,
    MAX_ETA,
)

import numpy as np
from selection.lumi import get_lumi_mask


class NanoReducer:

    def __init__(
        self,
        input_file,
        jet_selection=True,
        electron_selection=True,
        muon_selection=True,
        photon_selection=True,
        event_selection=True,
        apply_lumimask=True,
        lumimask_json=None,
        apply_trigger=True,
        trigger_paths=None,
    ):
        self.reader = NanoReader(input_file)
        self.jet_selection = jet_selection
        self.electron_selection = electron_selection
        self.muon_selection = muon_selection
        self.photon_selection = photon_selection
        self.event_selection = event_selection

        self.apply_lumimask = apply_lumimask
        self.lumimask_json = lumimask_json
        self.apply_trigger = apply_trigger
        self.trigger_paths = (
            trigger_paths
            if trigger_paths is not None
            else ["HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId"]
        )

    def run(self):

        collections = {}
        store = EventStore()

        # Detect MC vs Data
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

            sum_genw_orig = float(ak.sum(genWeight))
            n_events_orig = len(genWeight)

            store.add_metadata("sum_genw_presel", float(ak.sum(genWeight)))
            store.add_metadata("n_events_presel", len(genWeight))
        else:
            run_num = self.reader.read_scalar("run")
            original_index = ak.local_index(run_num)

            store.add_temp("original_index", original_index)
            store.add_metadata("n_events_presel", len(run_num))

        
        # Initial Cutflow Trackers
        n_all = len(original_index)
        store.add_metadata("cutflow_all", n_all)

        # 1: Lumi Mask Filtering (Data Only)
        # -------------------------------------------------------------
        if not is_mc and self.apply_lumimask and self.lumimask_json:
            run_arr = self.reader.read_scalar("run")
            lumi_arr = self.reader.read_scalar("luminosityBlock")
            lumi_mask = get_lumi_mask(
                ak.to_numpy(run_arr),
                ak.to_numpy(lumi_arr),
                self.lumimask_json,
            )
        else:
            lumi_mask = ak.ones_like(original_index, dtype=bool)

        store.add_metadata("cutflow_lumi_mask", int(ak.sum(lumi_mask)))

        # 2: HLT Trigger Filtering (Data Only)
        # -------------------------------------------------------------
        if not is_mc and self.apply_trigger:
            trigger_masks = []
            for trg in self.trigger_paths:
                try:
                    trg_val = self.reader.read_scalar(trg)
                    trigger_masks.append(trg_val)
                except Exception:
                    print(f"Trigger branch {trg} not found in input file.")

            if len(trigger_masks) > 0:
                hlt_mask = trigger_masks[0]
                for m in trigger_masks[1:]:
                    hlt_mask = hlt_mask | m
            else:
                hlt_mask = ak.ones_like(original_index, dtype=bool)
        else:
            hlt_mask = ak.ones_like(original_index, dtype=bool)
        # Combine Lumi + HLT
        lumi_hlt_mask = lumi_mask & hlt_mask
        store.add_metadata("cutflow_hlt", int(ak.sum(lumi_hlt_mask)))

        '''
        for collection in COLLECTIONS:
            print(f"Reading {collection}")
            collections[collection] = self.reader.read(collection)
        '''
        # -------------------------------------------------------------
        for collection in COLLECTIONS:
            print(f"Reading {collection}")
            obj = self.reader.read(collection)

            # keep defined fields, drop everything else
            if collection in KEEP_FIELDS:
                fields_to_keep = KEEP_FIELDS[collection]
                fields_to_drop = [f for f in obj.fields if f not in fields_to_keep]
                for field in fields_to_drop:
                    obj = ak.without_field(obj, field)

            collections[collection] = obj

        # 3: Photon Pre-mix selection
        photons = collections.get("Photon")

        if self.photon_selection and photons is not None:
            # At least 4 photons
            mask_4pho = lumi_hlt_mask & (ak.num(photons, axis=1) >= 4)
            store.add_metadata("cutflow_4photons", int(ak.sum(mask_4pho)))

            # Slicing first 4 photons
            pho4 = photons[:, :4]

            # Eta Acceptance
            abs_eta = np.abs(pho4.eta)
            eta_pass_each = (abs_eta < GAP_BARREL_ETA) | (
                (abs_eta > GAP_ENDCAP_ETA) & (abs_eta < MAX_ETA)
            )
            mask_eta = mask_4pho & ak.all(eta_pass_each, axis=1)
            store.add_metadata("cutflow_eta", int(ak.sum(mask_eta)))

            # PixelSeed Veto
            pixel_pass_each = pho4.pixelSeed == False
            mask_pixel = mask_eta & ak.all(pixel_pass_each, axis=1)
            store.add_metadata("cutflow_pixelSeed", int(ak.sum(mask_pixel)))

            # pT > 12 GeV
            pt_pass_each = pho4.pt > 12.0
            event_mask = mask_pixel & ak.all(pt_pass_each, axis=1)
            store.add_metadata("cutflow_pt_cuts", int(ak.sum(event_mask)))
        else:
            event_mask = lumi_hlt_mask

        # Save filtered branches to EventStore
        store.add_scalar("__original_index__", original_index[event_mask])

        for name, obj in collections.items():
            store.add_collection(name, obj[event_mask])

        for branch in SCALARS:
            print(f"Reading {branch}")
            try:
                value = self.reader.read_scalar(branch)
                store.add_scalar(branch, value[event_mask])
            except Exception:
                pass

        if is_mc:

            store.add_scalar(
                "sum_genw_presel", ak.full_like(genWeight[event_mask], sum_genw_orig)
            )

            for branch in SCALARS_MC:
                print(f"Reading MC scalar {branch}")
                try:
                    value = self.reader.read_scalar(branch)
                    store.add_scalar(branch, value[event_mask])
                except Exception:
                    pass
            
            for branch in WEIGHTS:
                print(f"Reading {branch}")
                try:
                    value = self.reader.read_weight(branch)
                    store.add_weight(branch, value[event_mask])
                except Exception:
                    pass

        return store
