import time
import os
import array
from ROOT import gROOT, TFile, TMVA, TCut, TTree
from h4g_tools.utils.runner_utils import get_parser_scan, handle_errors
from h4g_tools.utils.loading import pathSetup, loadSamples
from h4g_tools.utils.scales import Scales as ScalesCls

gROOT.SetBatch(True)

"""
Example usage:

python3 TMVA_test.py --gen -sig /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal_method1_SaS_Pileup_NFs_standardSelections_05Aug2025/ -bkg bdtIO/BDT_year_input_weights.parquet -y 2022postEE

"""



def gen_BDT(
    year: str,
    sig: str,
    bkg: str,
    debug: bool,
    xgb_name: str,
    eras: list = None,
    Scales: str = None,
    noRW: bool = False,
    mva_cut: float = None,
    EB_only: bool = False
) -> None:
    """
    Trains BDT over various masses for pseudoscalars using TMVA.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    print("\nLoading signal samples:")
    sig_samples = loadSamples(sig, eras=eras)
    if mva_cut is not None:
        assert type(mva_cut) is float, type(mva_cut)
        sig_samples = sig_samples[(sig_samples.pho3_mvaID >= mva_cut) & (sig_samples.pho4_mvaID >= mva_cut)]
    if EB_only:
        sig_samples = sig_samples[sig_samples.pho1_isScEtaEB & sig_samples.pho2_isScEtaEB & sig_samples.pho3_isScEtaEB & sig_samples.pho4_isScEtaEB]

    print("\nLoading background samples:")
    bkg_samples = loadSamples(bkg.replace("year", year), eras=eras)
    assert "Ndimreweight" in bkg_samples.columns
    if mva_cut is not None:
        assert type(mva_cut) is float, type(mva_cut)
        bkg_samples = bkg_samples[(bkg_samples.pho3_mvaID >= mva_cut) & (bkg_samples.pho4_mvaID >= mva_cut)]
    if EB_only:
        bkg_samples = bkg_samples[bkg_samples.pho1_isScEtaEB & bkg_samples.pho2_isScEtaEB & bkg_samples.pho3_isScEtaEB & bkg_samples.pho4_isScEtaEB]

    print("\nPreparing TTrees:")
    # Variable setup
    pho1_mvaID = array.array("d", [0])
    pho2_mvaID = array.array("d", [0])
    pho3_mvaID = array.array("d", [0])
    pho4_mvaID = array.array("d", [0])
    LeadPs_pt = array.array("d", [0])
    SubleadPs_pt = array.array("d", [0])
    dR_aa_mass_gggg = array.array("d", [0])
    LeadPs_interMass = array.array("d", [0])
    SubleadPs_interMass = array.array("d", [0])
    Ps_massDiff = array.array("d", [0])
    cos_ag = array.array("d", [0])
    weights = array.array("d", [0])

    # Signal branch setup
    sigTree = TTree()
    sigTree.Branch("pho1_mvaID", pho1_mvaID, "pho1_mvaID/D")
    sigTree.Branch("pho2_mvaID", pho2_mvaID, "pho2_mvaID/D")
    sigTree.Branch("pho3_mvaID", pho3_mvaID, "pho3_mvaID/D")
    sigTree.Branch("pho4_mvaID", pho4_mvaID, "pho4_mvaID/D")
    sigTree.Branch("LeadPs_pt", LeadPs_pt, "LeadPs_pt/D")
    sigTree.Branch("SubleadPs_pt", SubleadPs_pt, "SubleadPs_pt/D")
    sigTree.Branch("dR_aa_mass_gggg", dR_aa_mass_gggg, "dR_aa_mass_gggg/D")
    sigTree.Branch("LeadPs_interMass", LeadPs_interMass, "LeadPs_interMass/D")
    sigTree.Branch("SubleadPs_interMass", SubleadPs_interMass, "SubleadPs_interMass/D")
    sigTree.Branch("Ps_massDiff", Ps_massDiff, "Ps_massDiff/D")
    sigTree.Branch("cos_ag", cos_ag, "cos_ag/D")
    sigTree.Branch("weight", weights, "weight/D")

    # Fill signal TTree
    print("Filling signal tree...")
    for i in range(len(sig_samples)):
        sigTree.GetEntry(i)
        pho1_mvaID[0] = sig_samples["pho1_mvaID"][i]
        pho2_mvaID[0] = sig_samples["pho2_mvaID"][i]
        pho3_mvaID[0] = sig_samples["pho3_mvaID"][i]
        pho4_mvaID[0] = sig_samples["pho4_mvaID"][i]
        LeadPs_pt[0] = sig_samples["LeadPs_pt"][i]
        SubleadPs_pt[0] = sig_samples["SubleadPs_pt"][i]
        dR_aa_mass_gggg[0] = sig_samples["dR_aa_mass_gggg"][i]
        LeadPs_interMass[0] = sig_samples["LeadPs_interMass"][i]
        SubleadPs_interMass[0] = sig_samples["SubleadPs_interMass"][i]
        Ps_massDiff[0] = sig_samples["Ps_massDiff"][i]
        cos_ag[0] = sig_samples["cos_ag"][i]
        weights[0] = sig_samples["weight"][i]
        sigTree.Fill()

    # Background branch setup
    bkgTree = TTree()
    bkgTree.Branch("pho1_mvaID", pho1_mvaID, "pho1_mvaID/D")
    bkgTree.Branch("pho2_mvaID", pho2_mvaID, "pho2_mvaID/D")
    bkgTree.Branch("pho3_mvaID", pho3_mvaID, "pho3_mvaID/D")
    bkgTree.Branch("pho4_mvaID", pho4_mvaID, "pho4_mvaID/D")
    bkgTree.Branch("LeadPs_pt", LeadPs_pt, "LeadPs_pt/D")
    bkgTree.Branch("SubleadPs_pt", SubleadPs_pt, "SubleadPs_pt/D")
    bkgTree.Branch("dR_aa_mass_gggg", dR_aa_mass_gggg, "dR_aa_mass_gggg/D")
    bkgTree.Branch("LeadPs_interMass", LeadPs_interMass, "LeadPs_interMass/D")
    bkgTree.Branch("SubleadPs_interMass", SubleadPs_interMass, "SubleadPs_interMass/D")
    bkgTree.Branch("Ps_massDiff", Ps_massDiff, "Ps_massDiff/D")
    bkgTree.Branch("cos_ag", cos_ag, "cos_ag/D")
    bkgTree.Branch("weight", weights, "weight/D")

    # Fill background TTree
    print("Filling background tree...")
    for i in range(len(bkg_samples)):
        bkgTree.GetEntry(i)
        pho1_mvaID[0] = bkg_samples["pho1_mvaID"][i]
        pho2_mvaID[0] = bkg_samples["pho2_mvaID"][i]
        pho3_mvaID[0] = bkg_samples["pho3_mvaID"][i]
        pho4_mvaID[0] = bkg_samples["pho4_mvaID"][i]
        LeadPs_pt[0] = bkg_samples["LeadPs_pt"][i]
        SubleadPs_pt[0] = bkg_samples["SubleadPs_pt"][i]
        dR_aa_mass_gggg[0] = bkg_samples["dR_aa_mass_gggg"][i]
        LeadPs_interMass[0] = bkg_samples["LeadPs_interMass"][i]
        SubleadPs_interMass[0] = bkg_samples["SubleadPs_interMass"][i]
        Ps_massDiff[0] = bkg_samples["Ps_massDiff"][i]
        cos_ag[0] = bkg_samples["cos_ag"][i]
        if not noRW:
            weights[0] = bkg_samples["Ndimreweight"][i]
        else:
            weights[0] = 1.0

        bkgTree.Fill()

    if not os.path.exists(os.path.join(*(cwd, "outputs", "TMVA"))):
        os.mkdir(os.path.join(*(cwd, "outputs", "TMVA")))

    # TMVA setup
    TMVA.Tools.Instance()
    outputFile = TFile.Open(f"{cwd}/../models/tmva_classifier.root", "RECREATE")
    factory = TMVA.Factory(
        "4photon_classifier",
        outputFile,
        #"ROC=True, Silent=False, Color=True, AnalysisType='Classification'"
    )

    # Variables setup
    loader = TMVA.DataLoader("dataset")
    loader.AddVariable("pho1_mvaID")
    loader.AddVariable("pho2_mvaID")
    loader.AddVariable("pho3_mvaID")
    loader.AddVariable("pho4_mvaID")
    loader.AddVariable("LeadPs_pt")
    loader.AddVariable("SubleadPs_pt")
    loader.AddVariable("dR_aa_mass_gggg")
    loader.AddVariable("LeadPs_interMass")
    loader.AddVariable("SubleadPs_interMass")
    loader.AddVariable("Ps_massDiff")
    loader.AddVariable("cos_ag")

    signalWeight = bkgTree.GetEntries() / sigTree.GetEntries()
    backgroundWeight = 1.0

    loader.AddSignalTree(sigTree, signalWeight)
    loader.AddBackgroundTree(bkgTree, backgroundWeight)
    loader.SetSignalWeightExpression("weight")
    loader.SetBackgroundWeightExpression("weight")

    myCutS = TCut("")
    myCutB = TCut("")

    loader.PrepareTrainingAndTestTree(
        myCutS,
        myCutB,
        f"nTrain_Signal=0, nTrain_Background=0, SplitMode='Random', NormMode='NumEvents', V={debug}"
    )

    factory.BookMethod(
        loader,
        TMVA.Types.kBDT,
        "BDT",
        #"NTrees=30, MinNodeSize='2.5%', MaxDepth=3, BoostType='AdaBoost', AdaBoostBeta=0.5, UseBaggedBoost=True, BaggedSampleFraction=0.5, SeparationType='GiniIndex', nCuts=20"
        "NTrees=30, MinNodeSize='2.5%', BoostType=Grad, NegWeightTreatment=IgnoreNegWeightsInTraining, MaxDepth=3, UseRandomisedTrees=True, MinNodeSize=2, nCuts=20"
    )

    factory.TrainAllMethods()
    factory.TestAllMethods()
    factory.EvaluateAllMethods()

    outputFile.Close()


if __name__ == "__main__":
    parser = get_parser_scan()
    parser.add_argument("-noRW", "--noRW", default=False, action="store_true", help="Skip reweighting. Useful for testing efficacy in BDT training.")
    parser.add_argument("-mvaCut", "--mvaCut", default=-1.0, type=float, required=False, help="Cut to apply on photon MVA ID score.")
    parser.add_argument("-EB", "--EB_only", default=False, action="store_true", help="Use only events with all 4 photons in EB.")
    args = parser.parse_args()
    
    # Catch parser issues
    handle_errors(
        ["Need training samples.", args.sigInput is None, args.bkgInput is None, args.data is None],
        ["Must provide a year for processing.", args.year is None]
    )

    pathSetup(args.year)

    # Scan over masses 15-60 GeV in steps of 5 GeV
    masses = [x for x in range(15, 65, 5)]
    
    # Scales instance setup
    Scales = ScalesCls(args.year) 

    if args.year == "2022preEE":
        eras={
            "data": ["Run2022C", "Run2022D"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D"],
            "signal": [f"Signal_{m}_GeV_preEE" for m in masses]
        }
    elif args.year == "2022postEE":
        eras={
            "data": ["Run2022E", "Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_postEE" for m in masses]
        }
    elif args.year == "2022":
        eras={
            "data": ["Run2022C", "Run2022D", "Run2022E", "Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D", "RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_{prepost}" for m in masses for prepost in ["preEE", "postEE"]]
        }
    elif args.year == "2023preBPix":
        eras={
            "data": ["Run2023C"],
            "bkg": ["RunEvtMix2023C"],
            "signal": [f"Signal_{m}_GeV_preBPix" for m in masses]
        }
    elif args.year == "2022postBPix":
        eras={
            "data": ["Run2022D"],
            "bkg": ["RunEvtMix2022D"],
            "signal": [f"Signal_{m}_GeV_postBPix" for m in masses]
        }
    elif args.year == "2023":
        eras={
            "data": ["Run2023C", "Run2023D"],
            "bkg": ["RunEvtMix2023C", "RunEvtMix2023D"],
            "signal": [f"Signal_{m}_GeV_{prepost}" for m in masses for prepost in ["preBPix", "postBPix"]]
        }
    else:
        assert False, "Must specify year in Scales class!"

    start_time = time.time()
    gen_BDT(args.year, args.sigInput, args.bkgInput, args.debug, args.xgb_name, eras = eras, Scales=Scales, noRW=args.noRW, mva_cut=args.mvaCut, EB_only=args.EB_only)

    elapsed_time = time.time() - start_time
    print(f"Total elapsed time: {elapsed_time}s")




