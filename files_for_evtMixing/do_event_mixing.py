import os
import sys
import json
import ROOT
import argparse
import logging
import itertools
import numpy as np
from typing import ClassVar
from array import array as array
from pathlib import Path


# Make class for Tree creation here!
class TreeHelper:
    """
    Helper class to link TTree branches and variables during intermediate steps before Fill.
    """

    def __init__(
        self,
        tree: ROOT.TTree,
        access: ClassVar,
        branches: list,
        logger: logging.Logger,
        macro: str,
    ) -> None:
        """
        Variables are already set as data members of the TTree object, so no need to create anything here.
        """

        # List of all branches in NanoAOD v11: https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv11/2022postEE/doc_WZ_TuneCP5_13p6TeV_pythia8_Run3Summer22EENanoAODv11-126X_mcRun3_2022_realistic_postEE_v1-v1.html

        # Instance attributes containing arrays for branches used to be here. NOT needed since we can retain the TBranch structure from nanoAOD via CloneTree(0)
        self.tree = tree
        self.access = access  # Instance of Events class from MakeClass macro. Used to set branch values!
        ROOT.SetOwnership(self.tree, False)
        ROOT.SetOwnership(self.access, False)

        # Add photon branches so we can write our new values! All photon-specific branches (need size = 4)
        # Don't shuffle the Idx variables since they're needed for cross-ref in Coffea
        self.branches = set([
            "nPhoton",
            "Photon_seediEtaOriX",
            "Photon_cutBased",
            "Photon_electronVeto",
            "Photon_hasConversionTracks",
            "Photon_isScEtaEB",
            "Photon_isScEtaEE",
            "Photon_mvaID_WP80",
            "Photon_mvaID_WP90",
            "Photon_pixelSeed",
            "Photon_seedGain",
            "Photon_seediPhiOriY",
            "Photon_vidNestedWPBitmap",
            "Photon_ecalPFClusterIso",
            "Photon_electronIdx",
            "Photon_energyErr",
            "Photon_energyRaw",
            "Photon_esEffSigmaRR",
            "Photon_esEnergyOverRawE",
            "Photon_eta",
            "Photon_etaWidth",
            "Photon_haloTaggerMVAVal",
            "Photon_hcalPFClusterIso",
            "Photon_hoe",
            "Photon_hoe_Tower",
            "Photon_hoe_PUcorr",
            "Photon_jetIdx",
            "Photon_mvaID",
            "Photon_pfChargedIso",
            "Photon_pfChargedIsoPFPV",
            "Photon_pfChargedIsoWorstVtx",
            "Photon_pfPhoIso03",
            "Photon_pfRelIso03_all_quadratic",
            "Photon_pfRelIso03_chg_quadratic",
            "Photon_phi",
            "Photon_phiWidth",
            "Photon_pt",
            "Photon_r9",
            "Photon_s4",
            "Photon_sieie",
            "Photon_sieip",
            "Photon_sipip",
            "Photon_superclusterEta",
            "Photon_trkSumPtHollowConeDR03",
            "Photon_trkSumPtSolidConeDR04",
            "Photon_x_calo",
            "Photon_y_calo",
            "Photon_z_calo",
        ])

        # Need to make this manually, or the lists will be the same in memory
        self.nObjList = {
            #"nboostedTau": [None, None],
            #"nCorrT1METJet": [None, None],
            "nElectron": [None, None],
            #"nFatJet": [None, None],
            #"nFsrPhoton": [None, None],
            #"nIsoTrack": [None, None],
            "nJet": [None, None],
            #"nL1EG": [None, None],
            #"nL1EtSum": [None, None],
            #"nL1Jet": [None, None],
            #"nL1Mu": [None, None],
            #"nL1Tau": [None, None],
            #"nLowPtElectron": [None, None],
            #"nProton_multiRP": [None, None],
            #"nMuon": [None, None],
            #"nPPSLocalTrack": [None, None],
            #"nSoftActivityJet": [None, None],
            #"nProton_singleRP": [None, None],
            #"nSubJet": [None, None],
            #"nTauProd": [None, None],
            #"nTau": [None, None],
            #"nTrigObj": [None, None],
            #"nOtherPV": [None, None],
            #"nSV": [None, None],
        }

        self.skip_branches = set([
            "nGenPart",
            "GenJetAK8",
            "nGenProton",
            "nSubGenJetAK8",
            "nGenVisTau",
            "nLHEPdfWeight",
            "nLHEReweightingWeight",
            "nLHEScaleWeight",
            "nPSWeight",
            "nLHEPart",
            "nGenDressedLepton",
            "nGenIsolatedPhoton",
        ])

        # Useful for testing on MC. Otherwise, it's not necessary. Fine to leave it in!
        self.mc = False
        if len([b for b in branches if "genPart" in b]) > 0:
            logger.warning("Running over MC sample!")
            self.mc = True

            for b in ["Photon_genPartFlav", "Photon_genPartIdx"]:
                if b in branches:
                    self.branches.add(b)

        # Just in case I overspecify the branches I want... Exclude them below!
        # Note: Sometimes this if-statement behaves weirdly... It will be caught be the assert if it does
        if len(self.branches) > len(branches):
            logger.warning(f"Removed from self.branches: {', '.join([b for b in self.branches if b not in branches])}.")
            self.branches = set([b for b in self.branches if b in branches])
        assert sorted(self.branches) == sorted(branches), (len(self.branches), len(branches))

        ## Fill dictionares with original/shuffled photon branches for debugging
        self.mixed_values = {}
        self.original_values = {}
        for b in self.branches:
            # Use empty list so we don't have to specify array type early!
            self.mixed_values.update({b: [None] * 4})
            self.original_values.update({b: [None] * 4})

        # Notify if macro is the MC one BUT not using MC based on gen branches
        if not self.mc and "MC" in macro:
            logger.error("Incorrectly using MC macro but there are no Gen branches present! Data requires a goldenJSON!")


    @staticmethod
    def SkimTree(
        tree: ROOT.TTree,
        logger: logging.Logger,
        noSkim: bool = False,
        golden: str = None,
        isMC: bool = False
    ) -> ROOT.TTree:
        """
        Skims TTree by applying a 4 photons cut and checking the HLT bit ("just as good" alternative to messy preselections + corrections business).
        """

        if not noSkim:
            logger.info("Applying 4 photon cut and checking HLT bit on input tree...")
        else:
            logger.info("Specified noSkim: Applying 4 photon cut on input tree...")

        # 4 photons cut
        # Note: It is actually sufficienct to do "Photon_pt[3] > -1", but I include all 4 just for clarity
        n_photons = "Photon_pt[0] > -1 && Photon_pt[1] > -1 && Photon_pt[2] > -1 && Photon_pt[3] > -1"
        
        # Check for HLT bit so we can increase likelihood of good events later!
        hlt = "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId == 1"

        if not noSkim and golden is not None:
            return tree.CopyTree(f"{n_photons} && {hlt}")
        else:
            return tree.CopyTree(f"{n_photons}")


    @staticmethod
    def convert(
        vin: ROOT.TLeaf,
        idx: int = None,
        tp: type = None
    ) -> ROOT.Int_t | ROOT.Bool_t | ROOT.Short_t | ROOT.Long_t | ROOT.Float_t:
        """
        Need to convert to the proper native Python type - they're consistent but just stupid sometimes
        """

        # List of ROOT datatypes is here: https://root.cern.ch/doc/v634/RtypesCore_8h.html

        # Table for handling ROOT types
        type_table = {
            #"Char_t": ROOT.Char_t,
            #"UChar_t": ROOT.UChar_t,
            "Char_t": ROOT.Int_t,
            "UChar_t": ROOT.Int_t,
            "Bool_t": ROOT.Bool_t,
            "Short_t": ROOT.Short_t,
            "UShort_t": ROOT.UShort_t,
            "Int_t": ROOT.Int_t,
            "UInt_t": ROOT.UInt_t,
            "Float_t": ROOT.Float_t,
            "Long_t": ROOT.Long_t,
            "ULong64_t": int,
            #"ULong64_t": np.ulonglong,  # I have no idea why it prefers an int... ROOT.ULong64_t doesn't exist in pyROOT, it seems.
        }

        if tp is None:
            for key, val in type_table.items():
                # Expecting vin.GetTypeName() as a string - compare to key
                if vin.GetTypeName() == key:
                    if idx is not None:
                        return val(vin.GetValue(idx))
                    else:
                        return val(vin.GetValue())
        else:
            if idx is not None:
                return tp(vin.GetValue(idx))
            else:
                return tp(vin.GetValue())



    def fillOtherBranches(
        self,
        evt: int,
        to_shuffle: ROOT.TTree,
        logger: logging.Logger,
    ) -> None:
        """
        Sets self.tree branch values for any unmixed/unshuffled branches from to_shuffle. Filling occurs outside this function.
        """

        # GetEntry to ensure we're on the correct event
        to_shuffle.GetEntry(evt)

        # Get n<obj> variable to determine lengths of corresponding object arrays for this event. This exclude Photon_* branches
        for nObj in self.nObjList.keys():

            leaf = to_shuffle.GetLeaf(nObj)
            if not leaf:
                # Branch isn't present in skimmed ROOT file; skip silently
                continue

            try:
                self.nObjList[nObj][0] = int(to_shuffle.GetLeaf(nObj).GetValue())
            except ReferenceError:
                logger.warning(f"ReferenceError: Event {evt} for {nObj} in fillOtherBranches.")
                # Skipping only this branch!
                continue
            
            self.nObjList[nObj][1] = set([br.GetName() for br in to_shuffle.GetListOfBranches() if f"{nObj[1:]}_" in br.GetName() and "HLT_" not in br.GetName() and "L1_" not in br.GetName() and "DST_" not in br.GetName()])

            # The if-statements catch all jet branches together, so they need to be fixed before setting
            if nObj == "nFatJet":
                self.nObjList[nObj][1] = set([br for br in self.nObjList[nObj][1] if not br.startswith("Jet_") and not br.startswith("SubJet_") and not br.startswith("nL1Jet_")])
            elif nObj == "nSubJet":
                self.nObjList[nObj][1] = set([br for br in self.nObjList[nObj][1] if not br.startswith("Jet_") and not br.startswith("FatJet_") and not br.startswith("nL1Jet_")])
            elif nObj == "nJet":
                self.nObjList[nObj][1] = set([br for br in self.nObjList[nObj][1] if not br.startswith("FatJet_") and not br.startswith("SubJet_") and not br.startswith("nL1Jet_")])
            elif nObj == "nL1Jet":
                self.nObjList[nObj][1] = set([br for br in self.nObjList[nObj][1] if not br.startswith("Jet_") and not br.startswith("FatJet_") and not br.startswith("SubJet_")])

            # Same deal with tau branches
            if nObj == "nTau":
                self.nObjList[nObj][1] = set([br for br in self.nObjList[nObj][1] if not br.startswith("boostedTau_") and not br.startswith("L1Tau_")])
            elif nObj == "nboostedTau":
                self.nObjList[nObj][1] = set([br for br in self.nObjList[nObj][1] if not br.startswith("Tau_") and not br.startswith("L1Tau_")])
            elif nObj == "nL1Tau":
                self.nObjList[nObj][1] = set([br for br in self.nObjList[nObj][1] if not br.startswith("Tau_") and not br.startswith("boostedTau_")])
            logger.debug(f"{nObj}: {', '.join(self.nObjList[nObj][1])}.")
            
            # Write all branches coupled to n<obj> (of the form Obj_<var>)
            for b in self.nObjList[nObj][1]:
                logger.debug(f"Branch: {b}")

                # If the length of the array for this event is 0 then this loop will not execute - good behavior
                for index in range(self.nObjList[nObj][0]):
                    try:
                        getattr(self.access, b)[index] = TreeHelper.convert(to_shuffle.GetLeaf(b), index)
                    except IndexError as err:
                        logger.debug(err)

        # Now handle (uncoupled) branches of the form Obj_<var> (no index required). Omitting some HLT/L1 that are not of standard form and have multiple bits.
        # Uncoupled = list of branches per nObj (e.g. [[br1, br2, ...], ... [br25, ...]]) combined into one flat list like [br1, br2, ..., br25, ...]
        # Note: index 1 required for lambda since we put list of branches in that index of nObjList[nObj]
        uncoupled = set([br.GetName() for br in to_shuffle.GetListOfBranches() if br.GetName() not in list(itertools.chain(filter(None, map(lambda x: x[1], self.nObjList.values())))) and br.GetName() not in self.branches and br.GetName() not in self.skip_branches and "HLT" not in br.GetName() and "L1" not in br.GetName() and "DST" not in br.GetName()])
 
        # Add HLT/L1 paths here!
        for br in to_shuffle.GetListOfBranches():
            if "HLT_" in br.GetName() or "L1_" in br.GetName() or "Flag_" in br.GetName():
                uncoupled.add(br.GetName())

        # Set uncoupled branch values
        for b in uncoupled:
            # Note: Calling getattr(self.access, b) will sometimes yield an Overflow error for some reason, so don't do tests and print it!
            try:
                setattr(self.access, b, TreeHelper.convert(to_shuffle.GetLeaf(b))) 
            except TypeError as err:
                logger.debug(err)
            except RuntimeError as err:
                logger.debug(err)


    def shufflePhotons(
        self,
        event_to_write: int,
        event_to_write_cycle: int,
        event_to_read: int,
        to_shuffle: ROOT.TTree,
        index: int,
        logger: logging.Logger
    ) -> None:
        """
        Sets self.tree branch values to mixed/shuffled values from to_shuffle. Filling occurs outside this function.
        """

        # Get original photon values for all branches
        to_shuffle.GetEntry(event_to_write)
        for b in self.branches:
            self.original_values[b][index] = TreeHelper.convert(to_shuffle.GetLeaf(b), index)

        # Set event + offset for to_shuffle tree to GET mixed photons from to_shuffle
        to_shuffle.GetEntry(event_to_read, 1)

        # Get shuffled photon values for all branches
        for b in self.branches:
            if b != "nPhoton":
                self.mixed_values[b][index] = TreeHelper.convert(to_shuffle.GetLeaf(b), index)
            else:
                # Reset nPhoton branch to 4 since we have since skimmed the file!
                self.mixed_values[b] = ROOT.Int_t(4)
        
        # Set event for self.tree to SET shuffled photons in self.tree
        self.tree.GetEntry(event_to_write_cycle, 1)

        # Set value for photon index in an event (using mixed values conserves typing!)
        for key, val in self.mixed_values.items():
            if key != "nPhoton":
                getattr(self.access, key)[index] = self.mixed_values[key][index]
            else:
                # Manually convert to Int_t since nPhoton requires it. Always set this to 4 since we only keep the first 4 photons with mixing!
                setattr(self.access, key, ROOT.Int_t(4))


class SplitArgs(argparse.Action):
    """
    Splits comma-separated strings into lists for argparse.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(","))



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=str, action=SplitArgs, help="File to process with event mixing procedure. Can accept comma-separated paths for files.")
    parser.add_argument("outpath", type=str, help="Output path.")
    parser.add_argument("-g", "--goldenJSON", type=str, required=False, help="Golden JSON path. Only for data.")
    parser.add_argument("-nc", "--Ncycles", type=int, default=1, required=True, help="Event mixing cycles to run. Specifies offsets in steps of 3 to avoid overlap. Defaults to 1.")
    parser.add_argument("-offset", "--cycles_offset", type=int, default=0, required=False, help="Offset to use when running cycles for event mixing. Defaults to 0.")
    parser.add_argument("--do-3offset", action="store_true", help="Does offsets for cycles in sets of 3 instead of 1.")
    parser.add_argument("--no-event-mixing", action="store_true", help="Turns off event mixing procedure and only skims files. Useful for debugging or skimming data samples.")
    parser.add_argument("--debug", action="store_true", help="Run in debug mode. Limits to 20 events!")
    parser.add_argument("--show-all", action="store_true", help="Print all per-event information.")
    parser.add_argument("--noSkim", action="store_true", help="Do not skim samples with HLT prior to event mixing. Skims afterwards to reduce storage requirements.")
    args = parser.parse_args()

    # Want # of cycles to be at least one since we're just doing a loop over this number.
    assert args.Ncycles >= 1
    
    # Set this up here so the filled tree is only this size!
    if args.debug:
        debug_loop = 20

    # Use logger so log files are handled nicely
    logger = logging.getLogger()
    logging.basicConfig(stream=sys.stdout, format="%(asctime)s: %(levelname)s: %(message)s", level=logging.INFO if not args.debug else logging.DEBUG)
    logger.debug(" >>> Running in DEBUG MODE! <<<")
    logger.info(f"Processing file(s): {', '.join([os.path.abspath(f).split('/')[-1] for f in args.file])}")

    # Turn off RunRuntimeWarning from ROOT if not in debug mode
    if not args.debug:
        ROOT.gErrorIgnoreLevel = ROOT.kFatal
    else:
        logger.debug("All ROOT warnings are enabled due to debug mode.")

    # Load in relevant TTree
    logger.info("Loading input file(s)...")
    for idx, file in enumerate(args.file):
        if "root://" not in file:
            print(file)
            assert os.path.exists(file), os.path.abspath(file).split('/')[-1]
        else:
            logger.info(f"File #{idx} being loaded from remote using {file[file.find('root://')+7:file.find('///store')]}")

    # Use TChain so I don't have to handle a file directly! All branches are active at this point!
    # Also, so I can use multiple files for testing if I want.

    # Need to declare TChain in interpreter so I can pass it to the Events constructor (otherwise it tries to access a non-existent file -> null-pointer)
    ROOT.gInterpreter.Declare("""TChain* chain = new TChain("Events", "Events");""")
    raw_in_tree = ROOT.chain
    for file in args.file:
        retry = 3
        while retry > 0:
            raw_in_tree.AddFile(file)

            # Check if tree loaded properly
            if raw_in_tree.LoadTree(0) < 0 and retry != 1:
                logger.warning("File did not load properly! Retrying...")
                retry -= 1
            elif raw_in_tree.LoadTree(0) < 0 and retry == 1:
                logger.error("File did not load properly! Exiting!")
                raw_in_tree.Delete()
                sys.exit(1)
            elif raw_in_tree.LoadTree(0) >= 0:
                retry = 0

    if raw_in_tree.GetEntries() > 0:
        logger.info(f"Loaded input file with {raw_in_tree.GetEntries()} events.")
    else:
        logger.error(f"File did not load properly! Loaded input file with {raw_in_tree.GetEntries()} events. Exiting!")
        sys.exit(1)

    # Skip using MakeClass on input TTree
    # Overwriting in each job confuses ROOT, and we get an error.
    # Use compiled version instead so there are no issues with loading multiple times/simultaneously!
 
    # Only shuffle photon branches --> get names from raw_in_tree here!
    keepBranch = lambda k, *search: any([k.startswith(s) for s in search])
    branches = [key.GetName() for key in raw_in_tree.GetListOfBranches() if keepBranch(key.GetName(), "Photon_")]

    # Photon branches to exclude from shuffling so Coffea cross-ref works (wildcards * are used for context)
    to_exclude = []
    if len(to_exclude) > 0:
        logger.debug(f"Excluding any branches containing the following: {', '.join(to_exclude)}.")
    branches_ref = branches.copy()
    for b in branches:
        for exclude in to_exclude:
            if exclude.lower().replace("*","") in b.lower() and b in branches_ref:
                branches_ref.remove(b)
    branches_ref.append("nPhoton")
    branches = set(branches_ref)
    
    # Load in the macro. To create: 1. root <file>.root 2. Events->MakeClass() 3. Exit 4. root Events.C+
    # (not compiled b/c it's guaranteed that there are no syntax errors from MakeClass unless I introduce them!)
    cwd = os.path.dirname(os.path.abspath(__file__))
    macro = os.path.abspath(os.path.join(cwd, "Events_C.so"))  

    # Assume MC in this case! There is a confirmation/raised error in the TreeHelper constructor.
    if args.goldenJSON is None:
        macro = os.path.abspath(os.path.join(cwd, "EventsMC_C.so"))  
        
    assert os.path.exists(macro), macro

    # If I try to do this on the skimmed version, it loads "Memory Directory" since there is no storage location for this clone, then it fails!
    # Memory Directory came from doing MakeClass in this macro instead of prior to it. Still need to watch for errors, but that's where this comes from.
    # NOTE: Passing "chain->GetTree()" is critically important! Need this to load the skimmed file rather than whatever local tree Events_C.so derived from.
    # Also required to CloneTree(0) so we can write values to each event!
    ROOT.gROOT.LoadMacro(macro)
    ROOT.gInterpreter.Declare(f"Events cls(chain->GetTree()->CloneTree(0));")

    # Skim raw_in_tree so we know how many events to keep in filled_clone
    in_tree = TreeHelper.SkimTree(raw_in_tree, logger, args.noSkim, golden=args.goldenJSON)
    skimmed_events = in_tree.GetEntries()
    ROOT.SetOwnership(in_tree, False)
    logger.info(f"Skimmed down to {skimmed_events} events.")

    # Exit if skimmed event yield is zero since there's nothing more to do
    if skimmed_events == 0:
        logger.info(f"Finished processing {len(args.file)} file(s) with {skimmed_events} events! Exiting!")
        sys.exit(0)
        
    # Keep all branches - required to avoid multiple trees and saving issues later. Worth the I/O penalties
    ROOT.gInterpreter.Declare("TTree* empty_clone = cls.fChain->CloneTree(0);")

    # Set Events instance, empty clone, and filled clone as TreeHelper attributes.
    tree_helper = TreeHelper(ROOT.empty_clone, ROOT.cls, branches, logger, macro)
    assert tree_helper.tree.GetEntries() == 0, tree_helper.tree.GetEntries()

    # Don't include other trees/objects for consistency (they're not needed anyway).
    # If I do include them, it breaks for some reason, so just don't.

    logger.info("Preparing output tree...")
    # Get total # of events in skimmed input tree
    num_events = in_tree.GetEntries()
    # Inactive branches are not written to file so we keep all branches active!
    tree_helper.tree.SetBranchStatus("*", 1)

    # Start event loop here from 0 to num_events (exclusive). The full file is written so must shuffle the whole file!
    event_loop = num_events
    if args.debug:
        event_loop = debug_loop
        assert event_loop <= num_events, (event_loop, num_events)
    logger.info(f"Running event loop over {int(event_loop)} skimmed events:")
    
    # Warn if only skimming. Need to track this to ensure expected behavior
    if args.no_event_mixing:
        logger.warning("Specified NO event mixing. Filling new tree with skimmed events only!")

    logger.info(f"Entries in tree_helper before event loop: {tree_helper.tree.GetEntries()}")

    # Load Golden JSON
    golden_event_loop = []
    has_lumi_branches = hasattr(in_tree, "luminosityBlock") and hasattr(in_tree, "run") #addedByme
    if args.goldenJSON is not None and has_lumi_branches: #modifiedByme
        with open(args.goldenJSON, "r") as f:
            golden_json = json.load(f)

        for event in range(0, event_loop):
            # Check Golden JSON for Run and lumi section
            if str(in_tree.run) not in golden_json.keys():
                continue
            else:
                valid = False
                for subrange in golden_json[str(in_tree.run)]:
                    if in_tree.luminosityBlock in range(subrange[0], subrange[1]+1):
                        golden_event_loop.append(event)
                        valid = True
                if not valid:
                    continue
    else:
        golden_event_loop = [e for e in range(0, event_loop)]

    num_events = len(golden_event_loop)
    logger.info(f"{num_events} events remaining after Golden JSON check!")

    # Use "event_loop" for # of events to process and "num_events" for end of tree wrapping so a subset of events does not affect wrapping technique!
    skipped_cycles = 0
    for c in range(0 + args.cycles_offset, args.Ncycles + args.cycles_offset):
        if num_events == 0:
            continue

        logger.info(f" >> Processing cycle {c} <<")

        # Offsets/cycles are in steps of 1! Using steps of 3 makes the end of tree math much more complicated

        # Avoid additional cases by skipping the last 3 cycles before reaching length of event loop.
        # Otherwise, shuffled events will contain their original photons (e.g. event_loop - 3 has has photon 3 the same as original, and so on...) 
        cycle_check = None
        if args.do_3offset:
            # Steps of 3 for cycle offsets
            cycle_check = c > (golden_event_loop[-1] - 1)/3
        else:
            # Steps of 1 for cycle offsets
            cycle_check = c > golden_event_loop[-1] - 3 - 1

        if cycle_check:
            logger.warning(f"Proessing cycles past N events in tree will produce duplicate events! Skipping cycle {c}!")
            skipped_cycles += 1
            continue

        # Run over golden_event_loop
        #for event in range(0, event_loop):
        for idx, event in enumerate(golden_event_loop):
            # Each event was removed if it did not pass preselections (except diphotons) + 4 photons
            # Handle if loop is at the end of events. Set GetEntry(x, 1) to ensure read is only on active branches (i.e. all of them by design)!
            # Note: tree.GetEntry goes from 0 to tree.GetEntries() - 1

            # Set event for in_tree and tree_helper
            tree_helper.tree.GetEntry(event)
            in_tree.GetEntry(event)

            # Photon 1 is always the same. Need to offset write by cycle
            tree_helper.shufflePhotons(event, idx + c * event_loop, event, in_tree, 0, logger)           # Photon 1
            
            # Don't do event mixing if requested
            if args.no_event_mixing:
                tree_helper.shufflePhotons(event, idx + c * event_loop, event, in_tree, 1, logger)               # Photon 2
                tree_helper.shufflePhotons(event, idx + c * event_loop, event, in_tree, 2, logger)               # Photon 3
                tree_helper.shufflePhotons(event, idx + c * event_loop, event, in_tree, 3, logger)               # Photon 4

                # Set variables for unshuffled branches (use event NOT event_cycle since we're setting values from the input tree)
                logger.debug("Writing variables for unshuffled branches.")
                tree_helper.fillOtherBranches(event, in_tree, logger)

                # Skip the rest of this loop after filling tree
                tree_helper.tree.Fill()
                continue

            # Photon 4 end of range wrap.
            #if event == num_events - (2 + c) - 1:
            if idx == num_events - (2 + c) - 1:
                # Denote the event by the count from 1 rather than 0 to match the logger info
                if args.show_all:
                    logger.warning(f"Apply photon 4 end of range wrap (Event {num_events - (3 + c)}).")
                tree_helper.shufflePhotons(event, idx + c * event_loop, idx + 1 + c, in_tree, 1, logger)  # Photon 2
                tree_helper.shufflePhotons(event, idx + c * event_loop, idx + 2 + c, in_tree, 2, logger)  # Photon 3
                tree_helper.shufflePhotons(event, idx + c * event_loop, 0, in_tree, 3, logger)                  # Photon 4

            # Photon 3, 4 end of range wrap.
            #elif event == num_events - (1 + c) - 1:
            elif idx == num_events - (1 + c) - 1:
                # Denote the event by the count from 1 rather than 0 to match the logger info
                if args.show_all:
                    logger.warning(f"Apply photon 3,4 end of range wrap (Event {num_events - (2 + c)}).")
                tree_helper.shufflePhotons(event, idx + c * event_loop, idx + 1 + c, in_tree, 1, logger)  # Photon 2
                tree_helper.shufflePhotons(event, idx + c * event_loop, 0, in_tree, 2, logger)                  # Photon 3
                tree_helper.shufflePhotons(event, idx + c * event_loop, 1, in_tree, 3, logger)                  # Photon 4

            # Photon 2, 3, 4 end of range wrap.
            #elif event == num_events - (c) - 1:
            elif idx == num_events - (c) - 1:
                # Denote the event by the count from 1 rather than 0 to match the logger info
                if args.show_all:
                    logger.warning(f"Apply photon 2,3,4 end of range wrap (Event {num_events - (1 + c)}).")
                tree_helper.shufflePhotons(event, idx + c * event_loop, 0, in_tree, 1, logger)                  # Photon 2
                tree_helper.shufflePhotons(event, idx + c * event_loop, 1, in_tree, 2, logger)                  # Photon 3
                tree_helper.shufflePhotons(event, idx + c * event_loop, 2, in_tree, 3, logger)                  # Photon 4

            # Photon 2, 3, 4 end of range wrap for cycles.
            #elif event > num_events - (c) - 1:
            elif idx > num_events - (c) - 1:
                # Denote the event by the count from 1 rather than 0 to match the logger info
                if args.show_all:
                    logger.warning(f"Apply photon 2,3,4 end of range wrap with c (Event {event-1}).")
                tree_helper.shufflePhotons(event, idx + c * event_loop, event - (num_events - c - 1) + 0, in_tree, 1, logger)    # Photon 2
                tree_helper.shufflePhotons(event, idx + c * event_loop, event - (num_events - c - 1) + 1, in_tree, 2, logger)    # Photon 3
                tree_helper.shufflePhotons(event, idx + c * event_loop, event - (num_events - c - 1) + 2, in_tree, 3, logger)    # Photon 4

            # Otherwise, proceed as normal for bulk events.
            else:
                if args.show_all:
                    logger.debug("Normal mixing during bulk events.")
                tree_helper.shufflePhotons(event, idx + c * event_loop, idx + 1 + c, in_tree, 1, logger)  # Photon 2
                tree_helper.shufflePhotons(event, idx + c * event_loop, idx + 2 + c, in_tree, 2, logger)  # Photon 3
                tree_helper.shufflePhotons(event, idx + c * event_loop, idx + 3 + c, in_tree, 3, logger)  # Photon 4

            # Need to fill the events past the end of the original event_loop if Ncycles > 1
            tree_helper.tree.GetEntry(idx + c * event_loop)
            # Use tree_helper instance lists for convenience (even though they're the same as the values from the trees!)
            shuffled = ", ".join([f"{pt: .2f}".strip() for pt in list(tree_helper.mixed_values["Photon_pt"])])
            original = ", ".join([f"{pt: .2f}".strip() for pt in list(tree_helper.original_values["Photon_pt"])])
            shuffled = "[" + shuffled + "]"
            original = "[" + original + "]"

            # Keep this displayed so we can always confirm that event mixing is working!
            event_str = f"Event {event} / Cycle {c}"
            if args.show_all:
                logger.info(f"{event_str: <20} > Shuffled pt: {shuffled.strip():<30}\t Original pt: {original.strip():30}")
            assert list(tree_helper.mixed_values["Photon_pt"]) == list(tree_helper.tree.Photon_pt)[:4], (f"Event: {event}", list(tree_helper.mixed_values["Photon_pt"]), list(tree_helper.tree.Photon_pt)[:4])

            # Set variables for unshuffled branches (use event NOT event_cycle since we're setting values from the input tree)
            logger.debug("Writing variables for unshuffled branches.")
            tree_helper.fillOtherBranches(event, in_tree, logger)

            tree_helper.tree.Fill()


    print()
    if event_loop * (args.Ncycles - skipped_cycles) < in_tree.GetEntries() * (args.Ncycles - skipped_cycles):
        logger.warning(f"Processed fewer than total number of entries: {event_loop * (args.Ncycles - skipped_cycles)} < {in_tree.GetEntries() * (args.Ncycles - skipped_cycles)}")

    # Make directories as needed so ROOT can save outputs
    if "root://" not in args.outpath:
        rel_out = os.path.relpath(args.outpath, cwd)
        for i,p in enumerate(rel_out.split("/")[:-1]):
            # Ensure sub-directories are indeed SUB directories
            tmp = rel_out.split("/")[:i]
            tmp.append(p)

            if not os.path.exists(os.path.join(*tmp)):
                logger.debug(f"Creating directory {p} for outputs.")
                Path(os.path.join(*tmp)).mkdir(parents=True, exist_ok=True)
    else:
        # Directory will be handled by submit or resubmit scripts
        pass

    # Create empty out_file (at the end so there are no false positives for finished jobs)
    # Note the naming difference if skipping cycles!
    if int(tree_helper.tree.GetEntries()) > 0:
        true_outpath = args.outpath.replace(".root", f"_nc{args.Ncycles - skipped_cycles}_offset{args.cycles_offset}.root")
        logger.info(f"Saving outputs to file: {'/'.join(true_outpath.split('/')[-3:])}")
        out_file = ROOT.TFile.Open(true_outpath, "RECREATE")
        out_file.cd()
        tree_helper.tree.Write()

        out_file.Close()

        if "root://" not in true_outpath:
            assert os.path.exists(true_outpath), "Output writing failed!"
            assert os.stat(true_outpath).st_size != 0, f"Output file size is 0!"
    else:
        assert int(tree_helper.tree.GetEntries()) == 0
        logger.info("Not saving outputs to file since # of events is 0!")

    logger.info(f"Finished processing {len(args.file)} file(s) with {num_events} events! Wrote {int(tree_helper.tree.GetEntries())} events from {args.Ncycles - skipped_cycles} cycles.")
    if skipped_cycles > 0:
        logger.warning(f"Required to skip {skipped_cycles} cycles to avoid duplicate events!")

    print()
    # Access sample once per job! Otherwise, we may see xrootd access issues when trying to access the same file N times for all cycles in sequence
    # Output files contain all cycles!
