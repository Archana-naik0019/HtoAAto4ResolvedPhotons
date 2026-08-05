from tools.photonid_mva import (
    calculate_photonid_mva_run3,
    load_photonid_mva_run3,
)

from tools.SC_eta import add_photon_SC_eta
from tools.EELeak_region import veto_EEleak_flag
from tools.EcalBadCalibCrystal_events import remove_EcalBadCalibCrystal_events

from tools.flow_corrections import calculate_flow_corrections

from selections.photon_selections import (
    photon_preselection_h4g,
    addRhoCorrections,
)

from systematics import (
    object_systematics as available_object_systematics,
    object_corrections as available_object_corrections,
)

#----generic imports----
import functools
import operator
import os
import warnings
from typing import Any, Dict, List, Optional
import awkward
import numpy
import pandas
import sys
import vector
from coffea.nanoevents import NanoEventsFactory
from coffea import processor
from coffea.analysis_tools import Weights
from copy import deepcopy
import pyarrow
import pathlib
import shutil

import logging
#-----------------------
logger = logging.getLogger(__name__)
vector.register_awkward()
#------------------------


class HggBaseProcessor(processor.ProcessorABC):  # type: ignore
    def __init__(
        self,
        metaconditions: Dict[str, Any],
        systematics: Dict[str, List[str]] = None,
        corrections: Dict[str, List[str]] = None,
        year: Dict[str, List[str]] = None,
        doFlow_corrections: bool = False,
        validate_with_electrons: bool = False,
        output_format: str = "parquet",
    ) -> None:
        self.meta = metaconditions
        self.systematics = systematics if systematics is not None else {}
        self.corrections = corrections if corrections is not None else {}
        self.year = year if year is not None else {}
        self.doFlow_corrections = doFlow_corrections
        self.output_format = output_format
        self.validate_with_electrons = validate_with_electrons

        ###   Preselections   ###
        # eta
        self.gap_barrel_eta = 1.4442
        self.gap_endcap_eta = 1.566
        self.max_eta = 2.5

        # Electron veto
        #self.e_veto = 0.5
        self.e_veto = False

        ## HLT-mimicking selections
        # pT
        self.min_lead_pho_pt = 30.0
        self.min_sublead_pho_pt = 18.0

        # H/E
        self.hoe_eb = 0.08
        self.hoe_ee = 0.08

        # R9
        self.hlt_r9_eb = 0.5
        self.hlt_r9_ee = 0.8

        # Sigma_ieie
        #self.sigma_ieie_eb = 0.015
        #self.sigma_ieie_ee = 0.035
        self.sigma_ieie_eb = 0.011
        self.sigma_ieie_ee = 0.032

        # pfPhoIso + Rho corrections
        #self.pfPhoIso_eb = 4.0
        #self.pfPhoIso_ee = 4.0
        self.pfPhoIso_eb = 3.0
        self.pfPhoIso_ee = 3.0

        # TrackerIso
        self.trackerIso_eb = 6.0
        self.trackerIso_ee = 6.0

        # Multiple option isolation cuts
        self.iso_min_full5x5_r9 = 0.8
        self.iso_max_chad = 0.3
        self.iso_max_chad_rel = 20.0

        ###   Additional Preselections   ###
        self.diphoton_min_pt_lead = 30.0
        self.diphoton_min_pt_sublead = 18.0
        self.diphoton_min_lead_pt_mgg = 30.55 / 65  # ~0.47
        self.diphoton_min_sublead_pt_mgg = 18.20 / 65  # ~0.28

        # H4g-specific selections
        self.lead_pt = 30.0
        self.sublead_pt = 18.0
        self.min_pt = 15.0
        self.eta_gap_low = 1.4442
        self.eta_gap_high = 1.566
        self.eta_max = 2.5
        self.mass_range_low = 110.0
        self.mass_range_high = 180.0

        # EA values for Run3 from Egamma
        self.EA1_EB1 = 0.102056
        self.EA2_EB1 = -0.000398112
        self.EA1_EB2 = 0.0820317
        self.EA2_EB2 = -0.000286224
        self.EA1_EE1 = 0.0564915
        self.EA2_EE1 = -0.000248591
        self.EA1_EE2 = 0.0428606
        self.EA2_EE2 = -0.000171541
        self.EA1_EE3 = 0.0395282
        self.EA2_EE3 = -0.000121398
        self.EA1_EE4 = 0.0369761
        self.EA2_EE4 = -8.10369e-05
        self.EA1_EE5 = 0.0369417
        self.EA2_EE5 = -2.76885e-05

        logger.debug(f"Setting up processor with metaconditions: {self.meta}")

        #self.prefixes = {"pho_lead": "lead", "pho_sublead": "sublead"}
        self.prefixes = {
            "pho_lead": "lead",
            "pho_sublead": "sublead",
            "pho_subsublead": "subsublead",
            "pho_subsubsublead": "subsubsublead",
        }

        # initialize photonid_mva for Run3
        try:
            logger.info("Loading Run-3 PhotonID MVA models")

            photonid_dir = os.path.join(
                os.path.dirname(__file__),
                "tools",
                "flows",
            )

            self.photonid_mva = load_photonid_mva_run3(photonid_dir + "/")

        except Exception as e:
            warnings.warn(f"Could not instantiate Run-3 PhotonID MVA: {e}")
            self.photonid_mva = None

#-------------------------------------------------------------------
    def process(self, events: awkward.Array) -> Dict[Any, Any]:
    
            event_groups = {
                "two_photon_idx":   [4630, 4613, 4623],
                "three_photon_idx": [4615, 4617, 4607],
                "four_photon_idx":  [4609, 4636, 4608],
                "five_photon_idx":  [4606, 4611, 4655],
                "six_photon_idx":   [4786, 4826, 4876],
            }
    
            properties = [
                "pt",
                "eta",
                "r9",
                "pixelSeed",
                "hoe",
                "phi",
                "sieie",
                "trkSumPtHollowConeDR03",
                "trkSumPtSolidConeDR04",
                "pfPhoIso03",
            ]
    
            def output_write(tag, filename="event_dump.txt"):
                with open(filename, "a") as f:  # append to the same file
                    f.write("=" * 80 + "\n")
                    f.write(f"{tag}\n")
                    f.write("=" * 80 + "\n\n")
    
                    event_nums = events.event
    
                    for group_name, events_coll in event_groups.items():
                        f.write(group_name + "\n")
                        f.write("-" * len(group_name) + "\n\n")
    
                        for evt in events_coll:
                            idx = awkward.where(event_nums == evt)[0]
    
                            f.write(f"Event number : {evt}\n")
    
                            if len(idx) == 0:
                                f.write("Event not found.\n\n")
                                continue
    
                            idx = idx[0]
    
                            for prop in properties:
                                f.write(f"{prop:30s}: {events.Photon[prop][idx]*1}\n")
    
                            f.write("\n")
    
                        f.write("\n")
    
            dataset_name = events.metadata["dataset"]
            if "2018" in dataset_name:
                dataset_name = events.metadata["dataset"].replace("2018_", "")
    
            # Need Nevents for selection efficiency checks with signal samples
            # Number of initial events
            Nevents = {"initial_events": awkward.sum(events.genWeight) if "signal" in dataset_name.lower() else len(events)}
            Nphotons = {"initial_events": numpy.unique(awkward.num(events.Photon, axis=1).to_numpy(), return_counts=True)}
            Nphotons_EB = {"initial_events": numpy.unique(awkward.num(events.Photon[events.Photon.isScEtaEB], axis=1).to_numpy(), return_counts=True)}
            Nphotons_EE = {"initial_events": numpy.unique(awkward.num(events.Photon[events.Photon.isScEtaEE], axis=1).to_numpy(), return_counts=True)}
            Nphotons_OneEB = {"initial_events": numpy.unique(awkward.num(events.Photon[awkward.any(events.Photon.isScEtaEB, axis=1)], axis=1).to_numpy(), return_counts=True)}
    
            # data or monte carlo?
            self.data_kind = "mc" if hasattr(events, "GenPart") else "data"
    
            # here we start recording possible coffea accumulators
            # most likely histograms, could be counters, arrays, ...
            histos_etc = {}
            histos_etc[dataset_name] = {}
            if self.data_kind == "mc":
                histos_etc[dataset_name]["nTot"] = int(
                    awkward.num(events.genWeight, axis=0)
                )
                histos_etc[dataset_name]["nPos"] = int(awkward.sum(events.genWeight > 0))
                histos_etc[dataset_name]["nNeg"] = int(awkward.sum(events.genWeight < 0))
                histos_etc[dataset_name]["nEff"] = int(
                    histos_etc[dataset_name]["nPos"] - histos_etc[dataset_name]["nNeg"]
                )
                histos_etc[dataset_name]["genWeightSum"] = float(
                    awkward.sum(events.genWeight)
                )
            else:
                histos_etc[dataset_name]["nTot"] = int(len(events))
                histos_etc[dataset_name]["nPos"] = int(histos_etc[dataset_name]["nTot"])
                histos_etc[dataset_name]["nNeg"] = int(0)
                histos_etc[dataset_name]["nEff"] = int(histos_etc[dataset_name]["nTot"])
                histos_etc[dataset_name]["genWeightSum"] = float(len(events))
    
            if self.data_kind == "mc":
                # Add sum of gen weights before selection for normalisation in postprocessing
                histos_etc[dataset_name]["sum_genw_presel"] = str(awkward.sum(events.genWeight))
            else:
                histos_etc[dataset_name]["sum_genw_presel"] = "Data"
    
            # remove events affected by EcalBadCalibCrystal
            if self.data_kind == "data":
                excluded_years = ["2018", "2017", "2016preVFP", "2016postVFP"]
                if self.year[dataset_name][0] not in excluded_years:
                    events = remove_EcalBadCalibCrystal_events(events)
    
            # we need ScEta for corrections and systematics, it is present in NanoAODv13+ and can be calculated using PV for older versions
            events["Photon"] = add_photon_SC_eta(events.Photon, events.PV)
    
            # add veto EE leak branch for photons, could also be used for electrons
            if (
                self.year[dataset_name][0] == "2022EE"
                or self.year[dataset_name][0] == "2022postEE"
            ):
                events["Photon"] = veto_EEleak_flag(self, events.Photon)
    
            # read which systematics and corrections to process
            try:
                correction_names = self.corrections[dataset_name]
            except KeyError:
                correction_names = []
            try:
                systematic_names = self.systematics[dataset_name]
            except KeyError:
                systematic_names = []
    
            events["Photon"] = awkward.with_field(events.Photon, events.Photon.pt, "pt_raw")
            events["Electron"] = awkward.with_field(events.Electron, events.Electron.pt, "pt_raw")
    
            for correction_name in correction_names:
                if correction_name in available_object_corrections.keys():
                    logger.info(
                        f"Applying correction {correction_name} to dataset {dataset_name}"
                    )
                    varying_function = available_object_corrections[correction_name]
                    events = varying_function(
                        events=events, year=self.year[dataset_name][0]
                    )
                    continue
                else:
                    warnings.warn(f"Could not process correction {correction_name}.")
                    continue
    
            # Event mixing is done BEFORE processing with HiggsDNA!
            original_photons = events.Photon
    
            var_list = []  # so I get no errors for data
            # ***Computing the normalizing flow correction
            if self.data_kind == "mc" and self.doFlow_corrections:
    
                # Applying the Flow corrections to all photons before pre-selection
                counts = awkward.num(original_photons)
                corrected_inputs, var_list = calculate_flow_corrections(
                    original_photons,
                    events,
                    self.meta["flashggPhotons"]["flow_inputs"],
                    self.meta["flashggPhotons"]["Isolation_transform_order"],
                    year=self.year[dataset_name][0],
                )
    
                # Store the raw nanoAOD value and update photon ID MVA value for preselection
                original_photons["mvaID_nano"] = original_photons["mvaID"]
    
                for i in range(len(var_list)):
                    original_photons["raw_" + str(var_list[i])] = original_photons[str(var_list[i])]
                    original_photons[str(var_list[i])] = awkward.unflatten(corrected_inputs[:, i], counts)
    
                original_photons["mvaID"] = awkward.unflatten(
                    self.add_photonid_mva_run3(original_photons, events), counts
                )
    
            # systematic object variations
            for systematic_name in systematic_names:
                if systematic_name in available_object_systematics.keys():
                    systematic_dct = available_object_systematics[systematic_name]
                    if systematic_dct["object"] == "Photon":
                        logger.info(
                            f"Adding systematic {systematic_name} to photons collection of dataset {dataset_name}"
                        )
                        original_photons.add_systematic(
                            name=systematic_name,
                            kind=systematic_dct["args"]["kind"],
                            what=systematic_dct["args"]["what"],
                            varying_function=functools.partial(
                                systematic_dct["args"]["varying_function"],
                                events=events,
                                year=self.year[dataset_name][0],
                            ),
                        )
                    # to be implemented for other objects here
                else:
                    warnings.warn(
                        f"Could not process systematic variation {systematic_name}."
                    )
                    continue
    
            # Applying systematic variations
            photons_dct = {}
            photons_dct["nominal"] = original_photons
            logger.debug(original_photons.systematics.fields)
            for systematic in original_photons.systematics.fields:
                for variation in original_photons.systematics[systematic].fields:
                    # deepcopy to allow for independent calculations on photon variables with CQR
                    photons_dct[f"{systematic}_{variation}"] = deepcopy(
                        original_photons.systematics[systematic][variation]
                    )
    
            # Flatten and attach _down/_up
            variations_flattened = list(photons_dct.keys())
            variations_flattened.remove("nominal")
            variations = [item + suffix for item in variations_flattened for suffix in ["_down", "_up"]]
            variations.append("nominal")
            logger.debug(f"[systematics variations] {variations}")
    
            original_events = events
            for variation in variations:
                photons, jets = photons_dct["nominal"], events.Jet
                events = original_events
    
                if variation == "nominal":
                    pass
                elif variation in [*photons_dct]:
                    photons = photons_dct[variation]
    
                # recompute photonid_mva on the fly
                if self.photonid_mva:
                    counts = awkward.num(photons)
                    photons["mvaID"] = awkward.unflatten(self.add_photonid_mva_run3(photons, events), counts)
    
                # Apply rho correction to photons
                photons = addRhoCorrections(self, photons, events)
                assert "pfPhoIso03_rhoCorrected" in photons.fields
    
                # Minimum two photons needed to create diphotons
                photons = photons[awkward.num(photons, axis=1) >= 2]
                events = events[awkward.num(events.Photon, axis=1) >= 2]
    
                # Do pre-selections BEFORE making diphoton candidates
                photons, events, Nevents = photon_preselection_h4g(
                    self, photons, events, year=self.year[dataset_name][0], Nevents=Nevents
                )
                assert len(photons) == len(events), f"{len(photons)} {len(events)}"
    
                # Sort photons by pT
                photons = photons[awkward.argsort(photons.pt, axis=1, ascending=False)]
                assert len(photons) == len(events), f"{len(photons)} {len(events)}"
    
                # Keep only events with at least four selected photons
                mask = awkward.num(photons, axis=1) >= 4
                photons = photons[mask]
                events = events[mask]
    
                logger.info(f"Events with >=4 selected photons : {len(events)}")
    
                if len(events) == 0:
                    logger.warning("No events survive photon preselection.")
                    return histos_etc
    
                # Keep only the leading four photons
                photons = photons[:, :4]
    
                # output dataframe
                df = pandas.DataFrame({
                    "run": awkward.to_numpy(events.run),
                    "luminosityBlock": awkward.to_numpy(events.luminosityBlock),
                    "event": awkward.to_numpy(events.event),
                })
    
                # Save photon variables
                variables = list(dict.fromkeys([
                    "pt",
                    "eta",
                    "phi",
                    "r9",
                    "sieie",
                    "hoe",
                    "pfPhoIso03",
                    "pfPhoIso03_rhoCorrected",
                    "trkSumPtHollowConeDR03",
                    "trkSumPtSolidConeDR04",
                    "mvaID",
                    "mvaID_nano",
                ] + var_list))
    
                for i in range(4):
                    pho = photons[:, i]
                    for var in variables:
                        if var in pho.fields:
                            df[f"pho{i+1}_{var}"] = awkward.to_numpy(pho[var])
                        raw_name = f"raw_{var}"
                        if raw_name in pho.fields:
                            df[f"pho{i+1}_{raw_name}"] = awkward.to_numpy(pho[raw_name])
    
                # Save cutflow
                cutflow = (
                    pandas.Series(Nevents, name="events")
                    .rename_axis("selection")
                    .reset_index()
                )
    
                # Write output
                outdir = self.output_location if self.output_location else "."
                pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
                df.to_parquet(os.path.join(outdir, "photons.parquet"), index=False)
                cutflow.to_csv(os.path.join(outdir, "cutflow.csv"), index=False)
    
                logger.info(f"Saved {len(df)} events.")
                logger.info(f"Output directory : {outdir}")
    
                return histos_etc
#---------------------------------------------------------------------

    def postprocess(self, accumulant: Dict[Any, Any]) -> Any:
        pass

    def add_photonid_mva_run3(
        self, photons: awkward.Array, events: awkward.Array
    ) -> awkward.Array:

        photonid_mva_EB, photonid_mva_EE = self.photonid_mva

        rho = events.Rho.fixedGridRhoAll * awkward.ones_like(photons.pt)
        rho = awkward.flatten(rho)

        photons = awkward.flatten(photons)

        isEB = awkward.to_numpy(numpy.abs(photons.eta) < 1.5)
        mva_EB = calculate_photonid_mva_run3(
            [photonid_mva_EB, self.meta["flashggPhotons"]["inputs_EB"]], photons, rho
        )
        mva_EE = calculate_photonid_mva_run3(
            [photonid_mva_EE, self.meta["flashggPhotons"]["inputs_EE"]], photons, rho
        )
        mva = awkward.where(isEB, mva_EB, mva_EE)
        photons["mvaID_run3"] = mva

        return mva

    def apply_h4g_selection(
        self,
        photons: awkward.Array,
        events: awkward.Array,
        Nevents: dict,
    ):
        """
        Selections:
          - >= 4 selected photons
          - pT: 30/18/15/15 GeV
          - |eta| < 2.5
          - Remove EB-EE transition region
          - PixelSeed veto
          - 110 < m4g < 180 GeV
        """

        # Sort photons by pT
        photons = photons[awkward.argsort(photons.pt, axis=1, ascending=False)]

        # Require at least four photons
        mask = awkward.num(photons, axis=1) >= 4
        photons = photons[mask]
        events = events[mask]

        Nevents[">=4_photons"] = len(events)

        if len(events) == 0:
            return photons, events, Nevents

        # Keep only the leading four photons
        photons = photons[:, :4]

        # Four-photon pT cuts
        pt_mask = (
            (photons[:, 0].pt > self.lead_pt)
            & (photons[:, 1].pt > self.sublead_pt)
            & (photons[:, 2].pt > self.min_pt)
            & (photons[:, 3].pt > self.min_pt)
        )

        photons = photons[pt_mask]
        events = events[pt_mask]

        Nevents["pt_cuts_h4g"] = len(events)

        if len(events) == 0:
            return photons, events, Nevents

        # Eta acceptance
        eta_mask = (
            (abs(photons[:, 0].eta) < self.eta_max)
            & (abs(photons[:, 1].eta) < self.eta_max)
            & (abs(photons[:, 2].eta) < self.eta_max)
            & (abs(photons[:, 3].eta) < self.eta_max)
        )

        photons = photons[eta_mask]
        events = events[eta_mask]

        Nevents["eta_acceptance"] = len(events)

        if len(events) == 0:
            return photons, events, Nevents

        # EB-EE transition veto
        gap_mask = (
            ((abs(photons[:, 0].eta) < self.eta_gap_low) | (abs(photons[:, 0].eta) > self.eta_gap_high))
            & ((abs(photons[:, 1].eta) < self.eta_gap_low) | (abs(photons[:, 1].eta) > self.eta_gap_high))
            & ((abs(photons[:, 2].eta) < self.eta_gap_low) | (abs(photons[:, 2].eta) > self.eta_gap_high))
            & ((abs(photons[:, 3].eta) < self.eta_gap_low) | (abs(photons[:, 3].eta) > self.eta_gap_high))
        )

        photons = photons[gap_mask]
        events = events[gap_mask]

        Nevents["eta_gap_veto"] = len(events)

        if len(events) == 0:
            return photons, events, Nevents

        # Pixel seed veto
        pixel_mask = (
            (~photons[:, 0].pixelSeed)
            & (~photons[:, 1].pixelSeed)
            & (~photons[:, 2].pixelSeed)
            & (~photons[:, 3].pixelSeed)
        )

        photons = photons[pixel_mask]
        events = events[pixel_mask]

        Nevents["pixelSeed_veto"] = len(events)

        if len(events) == 0:
            return photons, events, Nevents

        # Four-photon invariant mass
        m4g = (
            photons[:, 0]
            + photons[:, 1]
            + photons[:, 2]
            + photons[:, 3]
        ).mass

        mass_mask = (
            (m4g > self.mass_range_low)
            & (m4g < self.mass_range_high)
        )

        photons = photons[mass_mask]
        events = events[mass_mask]

        Nevents["m4g_window"] = len(events)

        return photons, events, Nevents
