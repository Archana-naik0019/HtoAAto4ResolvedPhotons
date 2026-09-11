# flake8: noqa
# Suppresses flake8 on this file to pass pipeline

import xgboost
from higgs_dna.tools.chained_quantile import ChainedQuantileRegression
from higgs_dna.tools.diphoton_mva import calculate_diphoton_mva
from higgs_dna.tools.xgb_loader import load_bdt
from higgs_dna.tools.photonid_mva import calculate_photonid_mva, load_photonid_mva
from higgs_dna.tools.photonid_mva import calculate_photonid_mva_run3, load_photonid_mva_run3
from higgs_dna.tools.SC_eta import add_photon_SC_eta
from higgs_dna.tools.EELeak_region import veto_EEleak_flag
from higgs_dna.tools.EcalBadCalibCrystal_events import remove_EcalBadCalibCrystal_events
from higgs_dna.tools.gen_helpers import get_fiducial_flag, get_genJets, get_higgs_gen_attributes
from higgs_dna.tools.sigma_m_tools import compute_sigma_m
from higgs_dna.selections.photon_selections import photon_preselection_h4g, addRhoCorrections
from higgs_dna.selections.diphoton_selections import apply_fiducial_cut_det_level
from higgs_dna.selections.lepton_selections import select_electrons, select_muons
from higgs_dna.selections.jet_selections import select_jets, jetvetomap
from higgs_dna.selections.lumi_selections import select_lumis
from higgs_dna.utils.dumping_utils import (
    diphoton_ak_array,
    dump_ak_array,
    diphoton_list_to_pandas,
    ps_list_to_pandas,
    dump_pandas,
    get_obj_syst_dict,
)
from higgs_dna.utils.misc_utils import choose_jet
from higgs_dna.tools.flow_corrections import calculate_flow_corrections

from higgs_dna.tools.mass_decorrelator import decorrelate_mass_resolution

# from higgs_dna.utils.dumping_utils import diphoton_list_to_pandas, dump_pandas
from higgs_dna.metaconditions import photon_id_mva_weights
from higgs_dna.metaconditions import diphoton as diphoton_mva_dir
from higgs_dna.systematics import object_systematics as available_object_systematics
from higgs_dna.systematics import object_corrections as available_object_corrections
from higgs_dna.systematics import weight_systematics as available_weight_systematics
from higgs_dna.systematics import weight_corrections as available_weight_corrections

#from higgs_dna.utils.event_mixing import mix_photons_awkward # an attempt to integrate in evt-mixing after preselections

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
from coffea import processor
from coffea.analysis_tools import Weights
from copy import deepcopy
import pyarrow
import pathlib
import shutil

import logging

logger = logging.getLogger(__name__)

vector.register_awkward()


class HggBaseProcessor(processor.ProcessorABC):  # type: ignore
    def __init__(
        self,
        metaconditions: Dict[str, Any],
        systematics: Dict[str, List[str]] = None,
        corrections: Dict[str, List[str]] = None,
        apply_trigger: bool = True,
        output_location: str = None,
        taggers: List[Any] = None,
        nano_version: int = None,
        bTagEffFileName: str = None,
        trigger_group: str = ".*DoubleEG.*",
        analysis: str = "lowMassAnalysis",
        applyCQR: bool = False,
        skipJetVetoMap: bool = True,
        year: Dict[str, List[str]] = None,
        fiducialCuts: str = "classical",
        doDeco: bool = False,
        Smear_sigma_m: bool = False,
        doFlow_corrections: bool = False,
        validate_with_electrons: bool = False,
        output_format: str = "parquet",
    ) -> None:
        
        self.meta = metaconditions
        self.systematics = systematics if systematics is not None else {}
        self.corrections = corrections if corrections is not None else {}
        self.apply_trigger = apply_trigger
        self.nano_version = nano_version,
        self.bTagEffFileName = bTagEffFileName,
        self.output_location = output_location
        self.trigger_group = trigger_group
        self.analysis = analysis
        self.applyCQR = applyCQR,
        self.skipJetVetoMap = skipJetVetoMap
        self.year = year if year is not None else {}
        self.fiducialCuts = fiducialCuts
        self.doDeco = doDeco
        self.Smear_sigma_m = Smear_sigma_m
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
        self.hlt_r9_ee = 0.9

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

        self.taggers = []
        if taggers is not None:
            self.taggers = taggers
            self.taggers.sort(key=lambda x: x.priority)

        #self.prefixes = {"pho_lead": "lead", "pho_sublead": "sublead"}
        self.prefixes = {
            "pho_lead": "lead",
            "pho_sublead": "sublead",
            "pho_subsublead": "subsublead",
            "pho_subsubsublead": "subsubsublead",
        }
        if not self.doDeco:
            logger.info("Skipping Mass resolution decorrelation as required")
            #logger.info("Testing")
        else:
            logger.info("Performing Mass resolution decorrelation as required")

        # build the chained quantile regressions
        if numpy.array(self.applyCQR).flatten():
            try:
                self.chained_quantile: Optional[
                    ChainedQuantileRegression
                ] = ChainedQuantileRegression(**self.meta["PhoIdInputCorrections"])
            except Exception as e:
                warnings.warn(f"Could not instantiate ChainedQuantileRegression: {e}")
                self.chained_quantile = None
        else:
            logger.info("Skipping CQR as required")
            self.chained_quantile = None

        # initialize photonid_mva
        photon_id_mva_dir = os.path.dirname(photon_id_mva_weights.__file__)
        try:
            logger.debug(
                f"Looking for {self.meta['flashggPhotons']['photonIdMVAweightfile_EB']} in {photon_id_mva_dir}"
            )
            self.photonid_mva_EB = load_photonid_mva(
                os.path.join(
                    photon_id_mva_dir,
                    self.meta["flashggPhotons"]["photonIdMVAweightfile_EB"],
                )
            )
            self.photonid_mva_EE = load_photonid_mva(
                os.path.join(
                    photon_id_mva_dir,
                    self.meta["flashggPhotons"]["photonIdMVAweightfile_EE"],
                )
            )
        except Exception as e:
            warnings.warn(f"Could not instantiate PhotonID MVA on the fly: {e}")
            self.photonid_mva_EB = None
            self.photonid_mva_EE = None

        # initialize diphoton mva
        diphoton_weights_dir = os.path.dirname(diphoton_mva_dir.__file__)
        logger.debug(
            f"Base path to look for IDMVA weight files: {diphoton_weights_dir}"
        )

        try:
            self.diphoton_mva = load_bdt(
                os.path.join(
                    diphoton_weights_dir, self.meta["flashggDiPhotonMVA"]["weightFile"]
                )
            )
        except Exception as e:
            warnings.warn(f"Could not instantiate diphoton MVA: {e}")
            self.diphoton_mva = None

        # initialize h4g event selection BDT
        #h4g_bdt_dir = os.path.join(os.path.dirname(__file__), "/eos/user/a/arnaik/HiggsDNA_LCG/higgs-dna-2024/higgs_dna/metaconditions/h4g/")
        h4g_bdt_dir = "/eos/home-a/arnaik/HiggsDNA_LCG/higgs-dna-2024/higgs_dna/metaconditions/h4g/"

        try:
            if "2022" in str(year):
                bdt_year = 2022
            elif "2023" in str(year):
                bdt_year = 2023
            elif "2024" in str(year):
                bdt_year = 2024

            bdt = xgboost.XGBClassifier()
            self.h4g_bdt = load_bdt(
                os.path.join(
                    h4g_bdt_dir, f"ES_BDT_{bdt_year}.ubj"
                ),
                classifier=True
            )
            
            logger.info(
                f"Loaded h4g BDT: {h4g_bdt_dir}ES_BDT_{bdt_year}.ubj"
            )

            # !!!!! This is a required step! HiggsDNA can only load ubj, but my tools only output xgb... It's easy to convert, though.  !!!!
            #NOTE: To convert from xgb to ubj:
            logger.info("See comments in base.py on how to convert XGBoost model from .xgb to .ubj for use in HiggsDNA!")
            # Make sure you have an environment with XGBoost 2.1.1
            # Load the .xgb with bdt = xgboost.XGBClassifier() and bdt.load_model(<model>.xgb).
            # Save to .ubj with bdt.save_model(<model>.ubj) and put it in ../metaconditions/h4g/

        except Exception as e:
            warnings.warn(f"Could not instantiate h4g BDT: {e}")
            self.h4g_bdt = None

    def process_extra(self, meta: str, photons: awkward.Array, diphotons: awkward.Array, events: awkward.Array, signal: bool, variation: str, Nevents: bool) -> awkward.Array:
        # Run selections for pseudoscalars
        pseudos, diphotons, events, Nevents = self.produce_and_select_ps(
            meta,
            photons=photons,
            diphotons=diphotons,
            events=events,
            signal=signal,
            variation=variation,
            Nevents=Nevents
        )

        return pseudos, diphotons, events, Nevents

    def apply_filters_and_triggers(self, events: awkward.Array, Nevents: dict) -> awkward.Array: 
        if self.apply_trigger:
            trigger_names = []
            triggers = self.meta["TriggerPaths"][self.trigger_group][self.analysis]
            hlt = events.HLT
            for trigger in triggers:
                actual_trigger = trigger.replace("HLT_", "").replace("*", "")
                for field in hlt.fields:
                    if field == actual_trigger:
                        trigger_names.append(field)
            triggered = functools.reduce(
                operator.or_, (hlt[trigger_name] for trigger_name in trigger_names)
            )

        # Filters efficiency
        Nevents.update({"triggers": len(events[triggered]) if self.apply_trigger else len(events)})

        print(f"[{getattr(self, 'data_kind', 'unknown')}] Events before trigger: {len(events)} | Events after trigger: {len(events[triggered]) if self.apply_trigger else len(events)}")

        return events[triggered] if self.apply_trigger else events
    

                
    def process(self, events: awkward.Array) -> Dict[Any, Any]:
        
        #***********************************************************************************

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
            dataset_name = events.metadata["dataset"].replace("2018_","")

        # Need Nevents for selection efficiency checks with signal samples
        # Number of initial events
        Nevents = {"initial_events": awkward.sum(events.genWeight) if "signal" in dataset_name.lower() else len(events)}
        
        Nphotons = {"initial_events": numpy.unique(awkward.num(events.Photon, axis=1).to_numpy(), return_counts=True)}
        Nphotons_EB = {"initial_events": numpy.unique(awkward.num(events.Photon[events.Photon.isScEtaEB], axis=1).to_numpy(), return_counts=True)}
        Nphotons_EE = {"initial_events": numpy.unique(awkward.num(events.Photon[events.Photon.isScEtaEE], axis=1).to_numpy(), return_counts=True)}
        Nphotons_OneEB = {"initial_events": numpy.unique(awkward.num(events.Photon[awkward.any(events.Photon.isScEtaEB, axis=1)], axis=1).to_numpy(), return_counts=True)}

        # data or monte carlo? #temporarily changing this
        self.data_kind = "mc" if hasattr(events, "GenPart") or hasattr(events, "genWeight") else "data"

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

        # lumi mask
        if self.data_kind == "data":
            try:
                lumimask = select_lumis(self.year[dataset_name][0], events, logger)
                events = events[lumimask]
            except:
                logger.info(
                    f"[ lumimask ] Skip now! Unable to find year info of {dataset_name}"
                )

        # Lumi mask efficiency
        Nevents.update({"lumi_mask": awkward.sum(events.genWeight) if "signal" in dataset_name.lower() else len(events)})

        # No jets! skipJetVetoMap is set to True by default
        # apply jetvetomap: only retain events that without any jets in the EE leakage region
        if not self.skipJetVetoMap:
            events = jetvetomap(
                self, events, logger, dataset_name, year=self.year[dataset_name][0]
            )

        # metadata array to append to higgsdna output
        metadata = {}

        if self.data_kind == "mc":
            # Add sum of gen weights before selection for normalisation in postprocessing
            metadata["sum_genw_presel"] = str(awkward.sum(events.genWeight))
        else:
            metadata["sum_genw_presel"] = "Data"

        # Correct trigger has been applied so we can uncomment this AND the EcalBadCalibCrystal
        # apply filters and triggers

        #output_write("before trigger")
        
        if self.data_kind != "mc":
            self.apply_trigger = True
        else:
            self.apply_trigger = False

        events = self.apply_filters_and_triggers(events, Nevents)
        #output_write("after trigger, before corrections")
        # remove events affected by EcalBadCalibCrystal
        if self.data_kind == "data":
            excluded_years = ["2018", "2017", "2016preVFP", "2016postVFP"]
            if self.year[dataset_name][0] not in excluded_years:
                events = remove_EcalBadCalibCrystal_events(events)

        # we need ScEta for corrections and systematics, it is present in NanoAODv13+ and can be calculated using PV for older versions
        events["Photon"] = add_photon_SC_eta(events.Photon, events.PV)

        # Sometimes these can be saved to parquet as NaN entries. This can be dealt with in the to_pandas step.
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

        # If --Smear_sigma_m == True and no Smearing correction in .json for MC throws an error, since the pt scpectrum need to be smeared in order to properly calculate the smeared sigma_m_m
        if (
            self.data_kind == "mc"
            and self.Smear_sigma_m
            and "Smearing" not in correction_names
        ):
            warnings.warn(
                "Smearing should be specified in the corrections field in .json in order to smear the mass!"
            )
            sys.exit(0)

        # Since now we are applying Smearing term to the sigma_m_over_m i added this portion of code
        # specially for the estimation of smearing terms for the data events [data pt/energy] are not smeared!
        if self.data_kind == "data" and self.Smear_sigma_m:
            correction_name = "Smearing"

            logger.info(
                f"\nApplying correction {correction_name} to dataset {dataset_name}\n"
            )
            varying_function = available_object_corrections[correction_name]
            events = varying_function(events=events, year=self.year[dataset_name][0])

        events["Photon"] = awkward.with_field(events.Photon, events.Photon.pt, "pt_raw")
        events["Electron"] = awkward.with_field(events.Electron, events.Electron.pt, "pt_raw")
        
        for correction_name in correction_names:
            #print(correction_name, available_object_corrections.keys(), available_weight_corrections.keys())
            if correction_name in available_object_corrections.keys():
                logger.info(
                    f"Applying correction {correction_name} to dataset {dataset_name}"
                )
                varying_function = available_object_corrections[correction_name]
                events = varying_function(
                    events=events, year=self.year[dataset_name][0]
                )
            elif correction_name in available_weight_corrections:
                # event weight corrections will be applied after photon preselection / application of further taggers
                continue
            else:
                # may want to throw an error instead, needs to be discussed
                warnings.warn(f"Could not process correction {correction_name}.")
                continue
        #output_write("after corrections, before flow corrections")
        # Event mixing is done BEFORE processing with HiggsDNA! This procedure used to be done incorrectly here
        original_photons = events.Photon
        # NOTE: jet jerc systematics are added in the correction functions and handled later
        original_jets = events.Jet

        # Computing the normalizing flow correction
        if self.data_kind == "mc" and self.doFlow_corrections:

            # Applyting the Flow corrections to all photons before pre-selection
            counts = awkward.num(original_photons)
            corrected_inputs,var_list = calculate_flow_corrections(original_photons, events, self.meta["flashggPhotons"]["flow_inputs"], self.meta["flashggPhotons"]["Isolation_transform_order"], year=self.year[dataset_name][0])

            # Store the raw nanoAOD value and update photon ID MVA value for preselection
            original_photons["mvaID_nano"] = original_photons["mvaID"]

            # Store the raw values of the inputs and update the input values with the corrections since some variables used in the preselection
            for i in range(len(var_list)):
                original_photons["raw_" + str(var_list[i])] = original_photons[str(var_list[i])]
                original_photons[str(var_list[i])] = awkward.unflatten(corrected_inputs[:,i] , counts)

            original_photons["mvaID"] = awkward.unflatten(self.add_photonid_mva_run3(original_photons, events), counts)
        #output_write("after flow corrections, before processing start")
        # systematic object variations
        for systematic_name in systematic_names:
            if systematic_name in available_object_systematics.keys():
                systematic_dct = available_object_systematics[systematic_name]
                if systematic_dct["object"] == "Photon":
                    logger.info(
                        f"Adding systematic {systematic_name} to photons collection of dataset {dataset_name}"
                    )
                    original_photons.add_systematic(
                        # passing the arguments here explicitly since I want to pass the events to the varying function. If there is a more elegant / flexible way, just change it!
                        name=systematic_name,
                        kind=systematic_dct["args"]["kind"],
                        what=systematic_dct["args"]["what"],
                        varying_function=functools.partial(
                            systematic_dct["args"]["varying_function"],
                            events=events,
                            year=self.year[dataset_name][0],
                        )
                        # name=systematic_name, **systematic_dct["args"]
                    )
                # to be implemented for other objects here
            elif systematic_name in available_weight_systematics:
                # event weight systematics will be applied after photon preselection / application of further taggers
                continue
            else:
                # may want to throw an error instead, needs to be discussed
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
                ) #----I am perhaps encountering issue here while running the output from my event-mixing code, it works fine when systematics is left blank!

        # NOTE: jet jerc systematics are added in the corrections, now extract those variations and create the dictionary
        jerc_syst_list, jets_dct = get_obj_syst_dict(original_jets, ["pt", "mass"])
        # object systematics dictionary
        logger.debug(f"[ jerc systematics ] {jerc_syst_list}")

        # Build the flattened array of all possible variations
        variations_combined = []
        variations_combined.append(original_photons.systematics.fields)
        # NOTE: jet jerc systematics are not added with add_systematics
        variations_combined.append(jerc_syst_list)
        # Flatten
        variations_flattened = sum(variations_combined, [])  # Begin with empty list and keep concatenating
        # Attach _down and _up
        variations = [item + suffix for item in variations_flattened for suffix in ['_down', '_up']]
        # Add nominal to the list
        variations.append('nominal')
        logger.debug(f"[systematics variations] {variations}")

        original_events = events
        for variation in variations:
            photons, jets = photons_dct["nominal"], events.Jet
            events = original_events

            if variation == "nominal":
                pass  # Do nothing since we already get the unvaried, but nominally corrected objets above
            elif variation in [*photons_dct]:  # [*dict] gets the keys of the dict since Python >= 3.5
                photons = photons_dct[variation]
            elif variation in [*jets_dct]:
                jets = jets_dct[variation]
            do_variation = variation  # We can also simplify this a bit but for now it works

            if self.chained_quantile is not None:
                photons = self.chained_quantile.apply(photons, events)
            # recompute photonid_mva on the fly
            if self.photonid_mva_EB and self.photonid_mva_EE:
                # Issue with loading the non-run3 version b/c of # of inputs... Also not sure why it would even be the non-run3...
                #photons = self.add_photonid_mva(photons, events)
                # Updating mvaID for each variation (also with normalizing flows) in MC
                counts = awkward.num(photons)
                photons["mvaID"] = awkward.unflatten(self.add_photonid_mva_run3(photons, events), counts)

            # Apply rho correction to photons
            photons = addRhoCorrections(self, photons, events)
            assert "pfPhoIso03_rhoCorrected" in photons.fields 

            

            # Minimum two photons needed to create diphotons
            photons = photons[awkward.num(photons, axis=1) >= 2]
            events = events[awkward.num(events.Photon, axis=1) >= 2]  # apply to events for weight corrections later
            #output_write("just before preselections")
            
            

            
            # Do pre-selections BEFORE making diphoton candidates to remove bad events (on 1st/2nd photons only)
            photons = photons[awkward.argsort(photons.pt, axis=1, ascending=False)]
            
            photons, events, Nevents = photon_preselection_h4g(
                self, photons, events, year=self.year[dataset_name][0], Nevents=Nevents
            )
            assert len(photons) == len(events), f"{len(photons)} {len(events)}"

            
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            # Sort photons by pT (should already be by default unless event mixing!)
            # Event order remains unchanged which is why we don't have to sort events array

            #photons = photons[awkward.argsort(photons.pt, axis=1, ascending=False)]
            assert len(photons) == len(events), f"{len(photons)} {len(events)}"


            # Create diphotons and do selections before h4g-specfic selections
            #photons = photons[awkward.argsort(photons.pt, axis=1, ascending=False)]

            photons["event_number"] = events.event
            photons["run_number"] = events.run
            photons["luminosity_block"] = events.luminosityBlock
            photons["charge"] = awkward.zeros_like(
                photons.pt
            )  # added this because charge is not a property of photons in nanoAOD v11. We just assume every photon has charge zero...
            #diphotons = awkward.combinations(
            #    photons, 2, fields=["pho_lead", "pho_sublead"]
            #)
            #output_write("after preselections, before diphotons")
            diphotons = awkward.combinations(
                photons[:,:4], 4, fields=["pho_lead", "pho_sublead", "pho_subsublead", "pho_subsubsublead"]
            )
            
            
            
            #assert awkward.all(diphotons.pho_lead.pt) >= awkward.all(diphotons.pho_sublead.pt)
            assert (awkward.all(diphotons.pho_lead.pt) >= awkward.all(diphotons.pho_sublead.pt) >= awkward.all(diphotons.pho_subsublead.pt) >= awkward.all(diphotons.pho_subsubsublead.pt))
            # now turn the diphotons into candidates with four momenta and such
            diphoton_4mom = diphotons["pho_lead"] + diphotons["pho_sublead"]
            diphotons["pt"] = diphoton_4mom.pt
            diphotons["eta"] = diphoton_4mom.eta
            diphotons["phi"] = diphoton_4mom.phi
            diphotons["mass"] = diphoton_4mom.mass
            diphotons["charge"] = diphoton_4mom.charge
            diphotons["event_number"] = diphotons.pho_lead.event_number
            diphotons["run_number"] = diphotons.pho_lead.run_number
            diphotons["luminosity_block"] = diphotons.pho_lead.luminosity_block
            
            diphoton_pz = diphoton_4mom.z
            diphoton_e = diphoton_4mom.energy

            diphotons["rapidity"] = 0.5 * numpy.log(
                (diphoton_e + diphoton_pz) / (diphoton_e - diphoton_pz)
            )

            diphotons = awkward.with_name(diphotons, "PtEtaPhiMCandidate")

            diphotons = awkward.with_name(diphotons, "PtEtaPhiMCandidate")
            diphotons["pT1_m_gg"] = diphotons.pho_lead.pt / diphotons.mass 
            diphotons["pT2_m_gg"] = diphotons.pho_sublead.pt / diphotons.mass 
            if (self.data_kind == "mc") :
                # This is not used by FlashggFinalFits so we can set it to whatever we want, but it does need to exist!
                diphotons["dZ"] = events.GenVtx.z - events.PV.z if hasattr(events, "GenVtx") else awkward.ones_like(diphotons.pt)
                diphotons["weight"] = awkward.ones_like(diphotons.pt)
                diphotons["weight_central"] = awkward.ones_like(diphotons.pt)

            # Keep first (highest pT) diphotons - see AN 2016_410_v10, line 357
            diphotons = diphotons[awkward.argsort(diphotons.pt, axis=1, ascending=False)]
            #output_write("after diphotons argsort")
            
            diphotons = awkward.firsts(diphotons)
            dipho_mask = ~awkward.is_none(diphotons)
            assert len(dipho_mask) == len(photons), f"{len(dipho_mask)} {len(photons)}"
            
            # Ensure at least one diphoton candidate (this is per event)
            diphotons = diphotons[dipho_mask]
            photons = photons[dipho_mask]
            #output_write("after photons[dipho_mask]")
            events = events[dipho_mask]
            assert len(diphotons) == len(photons) == len(events), f"{len(diphotons)} {len(photons)} {len(events)}"

            # Pre-selections efficiency
            Nevents.update({"pre_selections": len(photons)})
            Nphotons.update({"pre_selections": numpy.unique(awkward.num(events.Photon, axis=1).to_numpy(), return_counts=True)})
            Nphotons_EB.update({"pre_selections": numpy.unique(awkward.num(events.Photon[events.Photon.isScEtaEB], axis=1).to_numpy(), return_counts=True)})
            Nphotons_EE.update({"pre_selections": numpy.unique(awkward.num(events.Photon[events.Photon.isScEtaEE], axis=1).to_numpy(), return_counts=True)})
            Nphotons_OneEB.update({"pre_selections": numpy.unique(awkward.num(events.Photon[awkward.any(events.Photon.isScEtaEB, axis=1)], axis=1).to_numpy(), return_counts=True)})

            # Pre-selections + 4 photon efficiency
            num_photons = awkward.num(photons, axis=1) >= 4
            photons = photons[num_photons]
            diphotons = diphotons[num_photons]
            events = events[num_photons]
            assert len(diphotons) == len(photons) == len(events), f"{len(diphotons)} {len(photons)} {len(events)}"
            Nevents.update({"pre_selections+4photon": len(photons)})
            Nphotons.update({"pre_selections+4photon": numpy.unique(awkward.num(events.Photon, axis=1).to_numpy(), return_counts=True)})
            Nphotons_EB.update({"pre_selections+4photon": numpy.unique(awkward.num(events.Photon[events.Photon.isScEtaEB], axis=1).to_numpy(), return_counts=True)})
            Nphotons_EE.update({"pre_selections+4photon": numpy.unique(awkward.num(events.Photon[events.Photon.isScEtaEE], axis=1).to_numpy(), return_counts=True)})
            Nphotons_OneEB.update({"pre_selections+4photon": numpy.unique(awkward.num(events.Photon[awkward.any(events.Photon.isScEtaEB, axis=1)], axis=1).to_numpy(), return_counts=True)})

            # --- NEW: dump state immediately after preselection, before h4g-specific selections ---
            if self.output_location is not None:
                self._dump_diphotons(diphotons, events, dataset_name, variation, subdir_tag="preselection")

            # H4g-specific processing
            #output_write("just before making pseudoscalars")
            
            pseudos, diphotons, events, Nevents = self.process_extra(events.metadata["dataset"], photons, diphotons, events, True if "signal" in dataset_name.lower() else False, variation, Nevents)
            
            assert len(pseudos) == len(diphotons) == len(events), f"{len(pseudos)} {len(diphotons)} {len(events)}"
            
            d = diphotons[diphotons["event_number"] == 4996]

            

            # Full selections efficiency
            Nevents.update({"selections": len(pseudos)})
            #output_write("just before BDT scoring")
            # Add H4g BDT scores to events
            if self.h4g_bdt is not None:
                def predictBDT(
                    bdt: xgboost.Booster,
                    ps: awkward.Array,
                ) -> awkward.Array:
                    features = pandas.DataFrame({
                        "pho1_mvaID": ps.pho1.mvaID.to_numpy().flatten(),
                        "pho2_mvaID": ps.pho2.mvaID.to_numpy().flatten(),
                        "pho3_mvaID": ps.pho3.mvaID.to_numpy().flatten(),
                        "pho4_mvaID": ps.pho4.mvaID.to_numpy().flatten(),
                        "LeadPs_pt": ps.LeadPs.pt.to_numpy().flatten(),
                        "SubleadPs_pt": ps.SubleadPs.pt.to_numpy().flatten(),
                        "dR_aa_mass_gggg": ps.dR_aa_mass_gggg.to_numpy().flatten(),
                        "LeadPs_interMass": ps.LeadPs_interMass.to_numpy().flatten(),
                        "SubleadPs_interMass": ps.SubleadPs_interMass.to_numpy().flatten(),
                        "Ps_massDiff": ps.Ps_massDiff.to_numpy().flatten(),
                        "cos_ag": ps.cos_ag.to_numpy().flatten(),
                    })

                    # Predictions are given as probability that input is of class "signal"
                    # Note: this is predict_proba b/c xgb_loader.py gives me a XGBClassifier since I hacked it.
        
                    results = bdt.predict_proba(features)
                    
                    # Add BDT scores to pseudoscalars array
                    ps["BDT_score"] = results[:,1]
                    logger.warning("BDT scores")
                    return ps

                pseudos = predictBDT(self.h4g_bdt, pseudos)

            # Apply first 4 photon cut to events since that's what is actually done in pseudos
            events["Photon"] = events.Photon[:,:4]
            Nphotons.update({"selections": numpy.unique(awkward.num(events.Photon, axis=1).to_numpy(), return_counts=True)})
            Nphotons_EB.update({"selections": numpy.unique(awkward.num(events.Photon[events.Photon.isScEtaEB], axis=1).to_numpy(), return_counts=True)})
            Nphotons_EE.update({"selections": numpy.unique(awkward.num(events.Photon[events.Photon.isScEtaEE], axis=1).to_numpy(), return_counts=True)})
            Nphotons_OneEB.update({"selections": numpy.unique(awkward.num(events.Photon[awkward.any(events.Photon.isScEtaEB, axis=1)], axis=1).to_numpy(), return_counts=True)})

            # Clean up and order Nevents for addition to metadata!
            #if "signal" in dataset_name.lower():
            Nevents_clean = {
                "initial_events": Nevents["initial_events"],
                "lumi_mask": Nevents["lumi_mask"],
                "triggers": Nevents["triggers"],
                "2photons": Nevents["2photons"],
                "pt_cuts": Nevents["pt_cuts"],
                "hoe": Nevents["hoe"],
                "r9": Nevents["r9"],
                "sigma_ieie": Nevents["sigma_ieie"],
                "photonIso": Nevents["photonIso"],
                "trackerIso": Nevents["trackerIso"],
                "eta": Nevents["eta"],
                "pixelSeed": Nevents["pixelSeed"],
                "miniAOD": Nevents["miniAOD"],
                "pre_selections": Nevents["pre_selections"],
                #"triggers": Nevents["triggers"],
                "pre_selections+4photon": Nevents["pre_selections+4photon"],
                "pt_cuts_pseudos": Nevents["pt_cuts_pseudos"],
                "eta_pseudos": Nevents["eta_pseudos"],
                "pixelSeed_pseudos": Nevents["pixelSeed_pseudos"],
                "mass": Nevents["mass"],
                "selections": Nevents["selections"],
            }

            diphotons["pho_lead"] = awkward.firsts(pseudos["pho1"])            
            diphotons["pho_sublead"] = awkward.firsts(pseudos["pho2"])
            diphotons["pho_subsublead"] = awkward.firsts(pseudos["pho3"])
            diphotons["pho_subsubsublead"] = awkward.firsts(pseudos["pho4"])
            
            

            # Omit BDT selections step so mass variations for bkg/data can be dealt with first!
            # Easy to apply later anyway...

            print()
            print()
            for k,v in Nevents_clean.items():
                print(k, v)
            print()
            print()
                
            # return if there is no surviving events (in this case, want at least one pseudoscalar)
            if len(pseudos) == 0:
                logger.info("No surviving events in this run, return now!")
                return histos_etc
            #output_write("before weight corrections")
            if self.data_kind == "mc":
                assert len(pseudos) == len(diphotons) == len(events), f"{len(pseudos)} {len(diphotons)} {len(events)}"
                # initiate Weight container here, after selection, since event selection cannot easily be applied to weight container afterwards
                event_weights = Weights(size=len(events),storeIndividual=True)
                # set weights to generator weights
                event_weights._weight = awkward.to_numpy(events["genWeight"])

                # corrections to event weights:
                for correction_name in correction_names:
                    if correction_name in available_weight_corrections:
                        logger.info(
                            f"Adding correction {correction_name} to weight collection of dataset {dataset_name}"
                        )
                        logger.warning(
                            f"KT -- Adding correction {correction_name} to weight collection of dataset {dataset_name}"
                        )
                        """
                        Some notes to ensure correct binning and corrections:
                        1. The diphotons array is used to get leading/subleading photon pt etc, e.g. PreselSF and TriggerSF
                        2. Otherwise, the events array is used to grab Photons and bin, e.g. scale/smearing or Pileup
                        3. The photonSFs use the events.Photon array (unlike the electronSFs)
                        """
                        common_args = {
                            "events": events,
                            "photons": diphotons,
                            "weights": event_weights,
                            "dataset_name": dataset_name,
                            "year": self.year[dataset_name][0],
                        }
                
                        varying_function = available_weight_corrections[correction_name]
                        event_weights = varying_function(**common_args)

                # systematic variations of event weights go to nominal output dataframe:
                if do_variation == "nominal":
                    for systematic_name in systematic_names:
                        if systematic_name in available_weight_systematics:
                            logger.info(
                                f"Adding systematic {systematic_name} to weight collection of dataset {dataset_name}"
                            )
                            logger.warning(
                                f"KT --- Adding systematic {systematic_name} to weight collection of dataset {dataset_name}"
                            )
                            if systematic_name == "LHEScale":
                                if hasattr(events, "LHEScaleWeight"):
                                    
                                    diphotons["nweight_LHEScale"] = awkward.num(
                                        events.LHEScaleWeight,
                                        axis=1,
                                    )
                                    diphotons[
                                        "weight_LHEScale"
                                    ] = events.LHEScaleWeight #adding original just in case need to revert, i wouldnt have to rerun
                                    

                                    # Renormalization scale variations, values based on LHEScaleweight docstring
                                    diphotons["weight_LHEScale_muR"] = events.LHEScaleWeight[:, [1, 4, 7]]
                                    diphotons["nweight_LHEScale_muR"] = awkward.num(
                                        events.LHEScaleWeight[:, [1, 4, 7]], axis=1
                                    )
                                    
                                    # Factorization scale variations
                                    diphotons["weight_LHEScale_muF"] = events.LHEScaleWeight[:, [3, 4, 5]]
                                    diphotons["nweight_LHEScale_muF"] = awkward.num(
                                        events.LHEScaleWeight[:, [3, 4, 5]], axis=1
                                    )
                                else:
                                    logger.info(
                                        f"No {systematic_name} Weights in dataset {dataset_name}"
                                    )
                            elif systematic_name == "LHEPdf":
                                if hasattr(events, "LHEPdfWeight"):
                                    # two AlphaS weights are removed
                                    diphotons["nweight_LHEPdf"] = (
                                        awkward.num(
                                            events.LHEPdfWeight,
                                            axis=1,
                                        )
                                        - 2
                                    )
                                    diphotons[
                                        "weight_LHEPdf"
                                    ] = events.LHEPdfWeight[
                                        :, :-2
                                    ]
                                else:
                                    logger.info(
                                        f"No {systematic_name} Weights in dataset {dataset_name}"
                                    )
                            else:
                                common_args = {
                                    "events": events,
                                    "photons": diphotons,
                                    "weights": event_weights,
                                    "dataset_name": dataset_name,
                                    "year": self.year[dataset_name][0],
                                }

                                varying_function = available_weight_systematics[systematic_name]
                                event_weights = varying_function(**common_args)

                diphotons["weight"] = event_weights.weight()
                diphotons["weight_central"] = event_weights.weight() / events["genWeight"]

                metadata["sum_weight_central"] = str(
                    awkward.sum(
                        event_weights.weight()
                    )
                )
                metadata["sum_weight_central_wo_bTagSF"] = str(
                    awkward.sum(
                        event_weights.weight() / event_weights.partial_weight(include=["bTagSF"])
                    )
                )

                # Handle variations
                if do_variation == "nominal":
                    if event_weights.variations:
                        logger.info(
                            "Adding systematic weight variations to nominal output file."
                        )
                    for modifier in event_weights.variations:
                        diphotons["weight_" + modifier] = event_weights.weight(modifier=modifier)

            # Add weight variables (=1) for data for consistent datasets
            else:
                diphotons["weight_central"] = awkward.ones_like(diphotons.pt)
                diphotons["weight"] = awkward.ones_like(diphotons.pt)

            # Add Nevents to metadata
            for key,val in Nevents_clean.items():
                # Need to convert to string, otherwise conversion to PyArrow Table fails for some reason
                metadata[key] = str(val)

            # Add Nphotons to metadata too
            for key, new_key in zip(Nphotons.keys(), ["Initial Events", "Pre-selections", "Pre-selections + At least 4 photons", "Pre-selections + At least 4 photons + Pseudoscalar selections"]):
                # Need to convert to string, otherwise conversion to PyArrow Table fails for some reason
                for N, weight in zip(Nphotons[key][0], Nphotons[key][1]):
                    metadata[f"Nphotons {new_key} {N}"] = str(weight)
                for N, weight in zip(Nphotons_EB[key][0], Nphotons_EB[key][1]):
                    metadata[f"Nphotons_EB {new_key} {N}"] = str(weight)
                for N, weight in zip(Nphotons_EE[key][0], Nphotons_EE[key][1]):
                    metadata[f"Nphotons_EE {new_key} {N}"] = str(weight)
                for N, weight in zip(Nphotons_OneEB[key][0], Nphotons_OneEB[key][1]):
                    metadata[f"Nphotons_OneEB {new_key} {N}"] = str(weight)

            if self.output_location is not None:
                assert ~awkward.any(awkward.is_none(diphotons)) and ~awkward.any(awkward.is_none(pseudos)), f"{~awkward.any(awkward.is_none(diphotons))} {~awkward.any(awkward.is_none(pseudos))}"
                assert len(diphotons) == len(pseudos)

                df = diphoton_list_to_pandas(self, diphotons)
                df_ps = ps_list_to_pandas(pseudos)

                final_df = pandas.concat([
                    df,
                    df_ps,
                    pandas.Series(events.PV.npvs.to_list(), name="npvs"),
                    pandas.Series(events.PV.npvsGood.to_list(), name="npvsGood"),
                ], axis=1)

                # Exclude awkward arrays inside the DataFrame. The ROOT merging script in flashgg FinalFits doesn't like that
                if self.data_kind == "mc":
                    final_df = pandas.concat([
                        final_df,
                        pandas.Series(events.HLT.Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId.to_list(), name="HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId"),
                        pandas.Series(events.genWeight.to_list(), name="genWeight"),
                        pandas.Series(events.Pileup.nPU.to_list(), name="nPU"),
                        pandas.Series(events.Pileup.nTrueInt.to_list(), name="nTrueInt"),
                    ], axis=1)

                # Drop EELeak columns
                to_remove = [col for col in final_df.columns if "EELeak" in col]
                final_df.drop(to_remove, inplace=True, errors="ignore")

                fname = (events.attrs["@events_factory"]._partition_key.replace("/", "_") + ".%s" % self.output_format)
                fname = (fname.replace("%2F","")).replace("%3B1","")
                subdirs = []
                if "dataset" in events.metadata:
                    subdirs.append(events.metadata["dataset"])
                subdirs.append(variation)

                final_arr = awkward.Array(pyarrow.Table.from_pandas(final_df))
                final_table = awkward.to_arrow_table(final_arr, extensionarray=False)

                # Add efficiency to metadata
                merged_metadata = {**metadata, **(final_table.schema.metadata or {})}
                final_table = final_table.replace_schema_metadata(merged_metadata)

                # Set up filename and (sub)directories
                local_file = (os.path.join(".", fname))
                merged_subdirs = os.path.sep.join(subdirs)
                destination = (os.path.join(self.output_location, os.path.join(merged_subdirs, fname)))

                # Write pyarrow table to parquet file for now
                pyarrow.parquet.write_table(final_table, local_file)
                logger.warning("Parquet saving started")
                logger.warning(f"Destination : {destination}")
                # Save table to parquet and move file to proper directories
                dirname = os.path.dirname(destination)
                if not os.path.exists(dirname):
                    pathlib.Path(dirname).mkdir(parents=True, exist_ok=True)
                shutil.copy(local_file, destination)
                assert os.path.isfile(destination)
                pathlib.Path(local_file).unlink()
                logger.warning("Parquet saving ended")

        return histos_etc

    def postprocess(self, accumulant: Dict[Any, Any]) -> Any:
        #raise NotImplementedError
        pass

    def add_diphoton_mva(
        self, diphotons: awkward.Array, events: awkward.Array
    ) -> awkward.Array:
        return calculate_diphoton_mva(
            (self.diphoton_mva, self.meta["flashggDiPhotonMVA"]["inputs"]),
            diphotons,
            events,
        )

    #----------------adding it to pool out files immediately after presel-----------
    def _dump_diphotons(self, diphotons, events, dataset_name, variation, subdir_tag):
        """
        Write out the current diphoton/event state to parquet.
        Reuses the same writing logic as the final output, but skips
        pseudoscalar-level info (which doesn't exist yet at preselection time).
        """
        df = diphoton_list_to_pandas(self, diphotons)
        final_df = pandas.concat([
            df,
            pandas.Series(events.PV.npvs.to_list(), name="npvs"),
            pandas.Series(events.PV.npvsGood.to_list(), name="npvsGood"),
        ], axis=1)
    
        if self.data_kind == "mc":
            final_df = pandas.concat([
                final_df,
                pandas.Series(events.genWeight.to_list(), name="genWeight"),
            ], axis=1)
    
        to_remove = [col for col in final_df.columns if "EELeak" in col]
        final_df.drop(to_remove, inplace=True, errors="ignore")
    
        fname = (events.attrs["@events_factory"]._partition_key.replace("/", "_") + ".%s" % self.output_format)
        fname = (fname.replace("%2F", "")).replace("%3B1", "")
    
        subdirs = []
        if "dataset" in events.metadata:
            subdirs.append(events.metadata["dataset"])
        subdirs.append(variation)
        subdirs.append(subdir_tag)  # keeps this dump from clashing with the final output
    
        final_arr = awkward.Array(pyarrow.Table.from_pandas(final_df))
        final_table = awkward.to_arrow_table(final_arr, extensionarray=False)
    
        local_file = os.path.join(".", fname)
        merged_subdirs = os.path.sep.join(subdirs)
        destination = os.path.join(self.output_location, os.path.join(merged_subdirs, fname))
    
        pyarrow.parquet.write_table(final_table, local_file)
        dirname = os.path.dirname(destination)
        if not os.path.exists(dirname):
            pathlib.Path(dirname).mkdir(parents=True, exist_ok=True)
        shutil.copy(local_file, destination)
        assert os.path.isfile(destination)
        pathlib.Path(local_file).unlink()
        logger.info(f"[{subdir_tag}] wrote {destination}")
    #---------------------------------------------------------

    def add_photonid_mva(
        self, photons: awkward.Array, events: awkward.Array
    ) -> awkward.Array:
        photons["fixedGridRhoAll"] = events.Rho.fixedGridRhoAll * awkward.ones_like(
            photons.pt
        )
        counts = awkward.num(photons, axis=-1)
        photons = awkward.flatten(photons)
        isEB = awkward.to_numpy(numpy.abs(photons.eta) < 1.5)
        mva_EB = calculate_photonid_mva(
            (self.photonid_mva_EB, self.meta["flashggPhotons"]["inputs_EB"]), photons
        )
        mva_EE = calculate_photonid_mva(
            (self.photonid_mva_EE, self.meta["flashggPhotons"]["inputs_EE"]), photons
        )
        mva = awkward.where(isEB, mva_EB, mva_EE)
        photons["mvaID"] = mva

        return awkward.unflatten(photons, counts)

    def add_photonid_mva_run3(
        self, photons: awkward.Array, events: awkward.Array
    ) -> awkward.Array:

        preliminary_path = os.path.join(os.path.dirname(__file__), '../tools/flows/run3_mvaID_models/')
        photonid_mva_EB, photonid_mva_EE = load_photonid_mva_run3(preliminary_path)

        rho = events.Rho.fixedGridRhoAll * awkward.ones_like(photons.pt)
        rho = awkward.flatten(rho)

        photons = awkward.flatten(photons)

        isEB = awkward.to_numpy(numpy.abs(photons.eta) < 1.5)
        mva_EB = calculate_photonid_mva_run3(
            [photonid_mva_EB, self.meta["flashggPhotons"]["inputs_EB"]], photons , rho
        )
        mva_EE = calculate_photonid_mva_run3(
            [photonid_mva_EE, self.meta["flashggPhotons"]["inputs_EE"]], photons, rho
        )
        mva = awkward.where(isEB, mva_EB, mva_EE)
        photons["mvaID_run3"] = mva

        return mva

    def add_corr_photonid_mva_run3(
        self, photons: awkward.Array, events: awkward.Array
    ) -> awkward.Array:

        preliminary_path = os.path.join(os.path.dirname(__file__), '../tools/flows/run3_mvaID_models/')
        photonid_mva_EB, photonid_mva_EE = load_photonid_mva_run3(preliminary_path)

        rho = events.Rho.fixedGridRhoAll * awkward.ones_like(photons.pt)
        rho = awkward.flatten(rho)

        photons = awkward.flatten(photons)

        # Now calculating the corrected mvaID
        isEB = awkward.to_numpy(numpy.abs(photons.eta) < 1.5)
        corr_mva_EB = calculate_photonid_mva_run3(
            [photonid_mva_EB, self.meta["flashggPhotons"]["inputs_EB_corr"]], photons, rho
        )
        corr_mva_EE = calculate_photonid_mva_run3(
            [photonid_mva_EE, self.meta["flashggPhotons"]["inputs_EE_corr"]], photons, rho
        )
        corr_mva = awkward.where(isEB, corr_mva_EB, corr_mva_EE)

        return corr_mva


    def produce_and_select_ps(self, meta: str, photons: awkward.Array, diphotons: awkward.Array, events: awkward.Array, signal: bool, variation: str, Nevents: dict) -> awkward.Array:
               
        # Sort photons by pt
        photons = photons[awkward.argsort(photons.pt, ascending=False, axis=1)]
        photons = photons[awkward.num(photons, axis=1) >= 4]

        # Set mass and charge fields to zeros
        photons["mass"] = numpy.zeros_like(photons.pt)
        photons["charge"] = numpy.zeros_like(photons.pt)
        
        # Make fields for photons
        photons = photons[awkward.num(photons, axis=1) >= 4]
        pseudos = awkward.combinations(photons[:,:4], 4, fields=["pho1", "pho2", "pho3", "pho4"])  # easiest/most compact method to add 4 fields
        pseudos["pho1"] = awkward.with_name(pseudos.pho1, "PtEtaPhiMCandidate")
        pseudos["pho2"] = awkward.with_name(pseudos.pho2, "PtEtaPhiMCandidate")
        pseudos["pho3"] = awkward.with_name(pseudos.pho3, "PtEtaPhiMCandidate")
        pseudos["pho4"] = awkward.with_name(pseudos.pho4, "PtEtaPhiMCandidate")
        pseudos["mass_gggg"] = (pseudos.pho1 + pseudos.pho2 + pseudos.pho3 + pseudos.pho4).mass

        # Add sumPt
        pseudos["sumPt_gggg"] = pseudos.pho1.pt + pseudos.pho2.pt + pseudos.pho3.pt + pseudos.pho4.pt

        # Pt cuts
        lead_pt = pseudos.pho1.pt > self.lead_pt
        sublead_pt = pseudos.pho2.pt > self.sublead_pt
        min_pt = (pseudos.pho3.pt > self.min_pt) & (pseudos.pho4.pt > self.min_pt)
        pseudos = pseudos[awkward.firsts(lead_pt & sublead_pt & min_pt)]
        diphotons = diphotons[awkward.firsts(lead_pt & sublead_pt & min_pt)]
        events = events[awkward.firsts(lead_pt & sublead_pt & min_pt)]
        assert len(pseudos) == len(diphotons) == len(events), (len(pseudos), len(diphotons), len(events))
        Nevents.update({"pt_cuts_pseudos": len(pseudos)})

        # Eta cuts
        eta_gap = ((abs(pseudos.pho1.eta) < self.eta_gap_low) | (abs(pseudos.pho1.eta) > self.eta_gap_high)) & ((abs(pseudos.pho2.eta) < self.eta_gap_low) | (abs(pseudos.pho2.eta) > self.eta_gap_high)) & ((abs(pseudos.pho3.eta) < self.eta_gap_low) | (abs(pseudos.pho3.eta) > self.eta_gap_high)) & ((abs(pseudos.pho4.eta) < self.eta_gap_low) | (abs(pseudos.pho4.eta) > self.eta_gap_high))
        eta_max = (abs(pseudos.pho1.eta) < self.eta_max) & (abs(pseudos.pho2.eta) < self.eta_max) & (abs(pseudos.pho3.eta) < self.eta_max) & (abs(pseudos.pho4.eta) < self.eta_max)
        eta_cuts = eta_gap & eta_max
        pseudos = pseudos[awkward.firsts(eta_cuts)]
        diphotons = diphotons[awkward.firsts(eta_cuts)]
        events = events[awkward.firsts(eta_cuts)]
        assert len(pseudos) == len(diphotons) == len(events), (len(pseudos), len(diphotons), len(events))
        Nevents.update({"eta_pseudos": len(pseudos)})

        # Pixel Seed cuts
        e_veto = (pseudos.pho1.pixelSeed == False) & (pseudos.pho2.pixelSeed == False) & (pseudos.pho3.pixelSeed == False) & (pseudos.pho4.pixelSeed == False)
        pseudos = pseudos[awkward.firsts(e_veto)]
        diphotons = diphotons[awkward.firsts(e_veto)]
        events = events[awkward.firsts(e_veto)]
        assert len(pseudos) == len(diphotons) == len(events), (len(pseudos), len(diphotons), len(events))
        Nevents.update({"pixelSeed_pseudos": len(pseudos)})

        # Mass cuts
        # NOTE: These cuts are per-event, while the previous are per-photon, hence only these ones need to be applied to diphotons and events. See assert checks above.
        mass_low = pseudos.mass_gggg > self.mass_range_low
        mass_high = pseudos.mass_gggg < self.mass_range_high
        #temporarily disbling mass cut
        mass_cuts = awkward.ones_like(awkward.any(mass_low & mass_high, axis=1), dtype=bool)
        #mass_cuts = awkward.any(mass_low & mass_high, axis=1)
        pseudos = pseudos[mass_cuts]
        Nevents.update({"mass": len(pseudos)})

        # Ensure no events without pseudoscalars (for all arrays)
        diphotons = diphotons[mass_cuts]
        events = events[mass_cuts]
        diphotons = diphotons[~awkward.is_none(pseudos)]
        events = events[~awkward.is_none(pseudos)]
        assert len(diphotons) == len(pseudos) == len(events), f"{len(diphotons)} {len(photons)} {len(events)}"

        assert len(pseudos) == len(pseudos[~awkward.is_none(pseudos)]), f"{len(pseudos)} {len(pseudos[~awkward.is_none(pseudos)])}. If there is None then at least one photon was removed. Removing doesn't trigger this!"
        psuedos = pseudos[~awkward.is_none(pseudos)]

        # Get length to compare after making pseudoscalars
        len_ps = len(pseudos)

        def makePsPermutation(ps_array, perm):
            phoA = ps_array[f"pho{perm[0]}"]
            phoB = ps_array[f"pho{perm[1]}"]
            assert awkward.all(awkward.all((phoA.pt >= phoB.pt), axis=1)), "PhoA pT must be bigger than PhoB pT!"

            # Ensure shape of array is (3, N) so we can keep the daughter photons AND the pseudoscalar together
            ps_perm = awkward.zip({"ps": phoA + phoB, "leading_pho": phoA, "subleading_pho": phoB})

            # NOTE: Old "lead_pho" method with ak.max yielded an nonsense 4-vector of the maximum elements from both photons...
            # This method here properly yields the leading/subleading decay photons from pseudoscalar candidates

            for field in ["ps", "leading_pho", "subleading_pho"]:
                ps_perm[field] = ps_perm[field]
                ps_perm[(field, "pt")] = ps_perm[field].pt
                ps_perm[(field, "eta")] = ps_perm[field].eta
                ps_perm[(field, "phi")] = ps_perm[field].phi
                ps_perm[(field, "mass")] = ps_perm[field].mass
                
            ps_perm["ps"] = awkward.with_name(ps_perm.ps, "PtEtaPhiMCandidate")
            ps_perm["leading_pho"] = awkward.with_name(ps_perm.leading_pho, "PtEtaPhiMCandidate")
            ps_perm["subleading_pho"] = awkward.with_name(ps_perm.subleading_pho, "PtEtaPhiMCandidate")

            return ps_perm

        # Start dM mixing procedure
        # Permutation A - (1,2,3,4)
        Ps11 = makePsPermutation(pseudos, (1,2))
        Ps12 = makePsPermutation(pseudos, (3,4))
        dM1 = abs(Ps11.ps.mass - Ps12.ps.mass)

        # Permutation B - (1,3,2,4)
        Ps21 = makePsPermutation(pseudos, (1,3))
        Ps22 = makePsPermutation(pseudos, (2,4))
        dM2 = abs(Ps21.ps.mass - Ps22.ps.mass)

        # Permutation C - (1,4,2,3)
        Ps31 = makePsPermutation(pseudos, (1,4))
        Ps32 = makePsPermutation(pseudos, (2,3))
        dM3 = abs(Ps31.ps.mass - Ps32.ps.mass)

        # Store eta information for daughter photons
        Ps11[("ps", "gg_inEB")] = (Ps11.leading_pho.eta < self.eta_gap_low) & (Ps11.subleading_pho.eta < self.eta_gap_low)
        Ps12[("ps", "gg_inEB")] = (Ps12.leading_pho.eta < self.eta_gap_low) & (Ps12.subleading_pho.eta < self.eta_gap_low)
        Ps21[("ps", "gg_inEB")] = (Ps21.leading_pho.eta < self.eta_gap_low) & (Ps21.subleading_pho.eta < self.eta_gap_low)
        Ps22[("ps", "gg_inEB")] = (Ps22.leading_pho.eta < self.eta_gap_low) & (Ps22.subleading_pho.eta < self.eta_gap_low)
        Ps31[("ps", "gg_inEB")] = (Ps31.leading_pho.eta < self.eta_gap_low) & (Ps31.subleading_pho.eta < self.eta_gap_low)
        Ps32[("ps", "gg_inEB")] = (Ps32.leading_pho.eta < self.eta_gap_low) & (Ps32.subleading_pho.eta < self.eta_gap_low)

        # Masks to determine lead/sublead pseudoscalars
        mask11 = Ps11.ps.pt >= Ps12.ps.pt
        mask12 = Ps11.ps.pt < Ps12.ps.pt
        mask21 = Ps21.ps.pt >= Ps22.ps.pt
        mask22 = Ps21.ps.pt < Ps22.ps.pt
        mask31 = Ps31.ps.pt >= Ps32.ps.pt
        mask32 = Ps31.ps.pt < Ps32.ps.pt

        # Minimize dM and create associated masks
        dM_min = awkward.min(awkward.concatenate((dM1, dM2, dM3), axis=1), axis=1)
        dM_mask1 = (dM_min == dM1)
        dM_mask2 = (dM_min == dM2)
        dM_mask3 = (dM_min == dM3)

        # Lead/sublead pseudoscalars for each permutation
        LeadPs1 = awkward.concatenate((Ps11[mask11], Ps12[mask12]), axis=1)
        SubleadPs1 = awkward.concatenate((Ps11[mask12], Ps12[mask11]), axis=1)
        LeadPs2 = awkward.concatenate((Ps21[mask21], Ps22[mask22]), axis=1)
        SubleadPs2 = awkward.concatenate((Ps21[mask22], Ps22[mask21]), axis=1)
        LeadPs3 = awkward.concatenate((Ps31[mask31], Ps32[mask32]), axis=1)
        SubleadPs3 = awkward.concatenate((Ps31[mask32], Ps32[mask31]), axis=1)

        # Apply masks and make final arrays
        dM = awkward.concatenate((dM1[dM_mask1], dM2[dM_mask2], dM3[dM_mask3]), axis=1)
        LeadPs = awkward.concatenate((LeadPs1[dM_mask1], LeadPs2[dM_mask2], LeadPs3[dM_mask3]), axis=1)
        SubleadPs = awkward.concatenate((SubleadPs1[dM_mask1], SubleadPs2[dM_mask2], SubleadPs3[dM_mask3]), axis=1)
        """
        # Testing the naive (12,34) mixing procedure
        dM = dM1
        LeadPs = LeadPs1
        SubleadPs = SubleadPs1
        """

        # Add dM, LeadPs, SubleadPs to events
        pseudos["dM"] = dM
        pseudos["LeadPs"] = LeadPs.ps
        pseudos[("LeadPs", "leading_pho")] = LeadPs.leading_pho
        pseudos[("LeadPs", "subleading_pho")] = LeadPs.subleading_pho
        pseudos["SubleadPs"] = SubleadPs.ps
        pseudos[("SubleadPs", "leading_pho")] = SubleadPs.leading_pho
        pseudos[("SubleadPs", "subleading_pho")] = SubleadPs.subleading_pho

        # Ensure length of pseudos stayed the same during mixing/optimization
        assert len_ps == len(pseudos)

        # Ensure fields in pseudos
        for field in ["pho1", "pho2", "pho3", "pho4", "LeadPs", "SubleadPs"]:
            pseudos[field] = pseudos[field]
            pseudos[(field, "pt")] = pseudos[field].pt
            pseudos[(field, "eta")] = pseudos[field].eta
            pseudos[(field, "phi")] = pseudos[field].phi
            pseudos[(field, "mass")] = pseudos[field].mass

            # Add iso_charged_rel to photons for plotting checks of pre-selections
            if field != "LeadPs" and field != "SubleadPs":
                pseudos[(field, "pfRelIso_photon_pt")] = (pseudos[field].pfRelIso03_chg if hasattr(pseudos[field], "pfRelIso03_chg") else pseudos[field].pfRelIso03_chg_quadratic) * pseudos[field].pt

        # Set behavior for pseudoscalar vectors (re vector module)
        pseudos["LeadPs"] = awkward.with_name(pseudos.LeadPs, "PtEtaPhiMCandidate")
        pseudos["SubleadPs"] = awkward.with_name(pseudos.SubleadPs, "PtEtaPhiMCandidate")
        pseudos["pT1_ma1"] = pseudos.LeadPs.leading_pho.pt / pseudos.LeadPs.mass
        pseudos["pT2_ma1"] = pseudos.LeadPs.subleading_pho.pt / pseudos.LeadPs.mass
        pseudos["pT1_ma2"] = pseudos.SubleadPs.leading_pho.pt / pseudos.SubleadPs.mass
        pseudos["pT2_ma2"] = pseudos.SubleadPs.subleading_pho.pt / pseudos.SubleadPs.mass

        # Store dR between pseudoscalars
        pseudos["dR_aa"] = pseudos.LeadPs.delta_r(pseudos.SubleadPs)
        pseudos["dR_aa_mass_gggg"] = pseudos.dR_aa / pseudos.mass_gggg
        
        # Store dR between daughter photons
        pseudos[("LeadPs", "dR_gg")] = pseudos.LeadPs.leading_pho.delta_r(pseudos.LeadPs.subleading_pho)
        pseudos[("SubleadPs", "dR_gg")] = pseudos.SubleadPs.leading_pho.delta_r(pseudos.SubleadPs.subleading_pho)

        # Make cos(theta_ag) variables
        # Get leading pseudoscalar in the lab frame such that the 4-vector is non-zero.
        # The negative sign is required by the code to undo the boost from the decaypseudos, diphotons, events, Nevents = self.process_extra(events.metadata["dataset"]
        pho1_Aframe = pseudos.LeadPs.leading_pho.boost(-pseudos.LeadPs.boostvec)  # Undoes the boost of the photon from the pseudoscalar
        pseudos["cos_ag"] = numpy.cos(pho1_Aframe.theta)
        # See slides from May 7, 2025 for confirmation on this behavior. This is the same as ROOT's implementation

        # Make hypMass flat distribution of nominal mass points for data/eventMixing
        if not signal:
            m_hyp = numpy.random.choice([float(x) for x in range(15, 65, 5)], len(pseudos.LeadPs.mass))
            assert len(m_hyp) == len(pseudos.LeadPs.mass)
        else:
            #m_hyp = numpy.full_like(pseudos.LeadPs.mass, float(meta.replace("2018_","")[7:9])) # using the three lines below instead
            import re
            mass_val = float(re.search(r"(\d+)", meta).group(1))
            m_hyp = numpy.full_like(pseudos.LeadPs.mass, mass_val)
            
            assert len(m_hyp) == len(pseudos.LeadPs.mass)

        # Make other BDT input variables
        pseudos["m_hyp"] = m_hyp
        #logger.warning(f"AT PHOTON ID SCALEFACTORS --- MASSPOINT : {m_hyp}")
        diphotons["masspoint"] = m_hyp #for custom PhotonID_SF values
        pseudos["LeadPs_interMass"] = (pseudos.LeadPs.mass - m_hyp) / pseudos.mass_gggg
        pseudos["SubleadPs_interMass"] = (pseudos.SubleadPs.mass - m_hyp) / pseudos.mass_gggg
        pseudos["Ps_massDiff"] = pseudos.LeadPs.mass - pseudos.SubleadPs.mass

        return pseudos, diphotons, events, Nevents

