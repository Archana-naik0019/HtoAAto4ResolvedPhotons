MAX_EVENTS_PER_FILE = 500000

COLLECTIONS = [
    #"Jet", #We won't need the whole collection in our analysis, but HiggsDNS crashes in its absence since the current base.py has jerc corrections defined (but never used). So may remove the jets collections from here if the 'jerc' correction lines are commented out from base.py
    "Photon",
    "Electron", #Not directly used in the analysis, but HiggsDNA crashes since the base.py stores electron raw pt (They are not used in the downstream analysis)
    "PuppiMET", # won't need the whole collection in our analysis
#    "PFMET",
    "PV",
#    "GenPart",
]

KEEP_FIELDS = {
    "Jet": [
        "pt",
        "eta",
        "phi",
        "neEmEF",
        "chEmEF", # these specific fields from the Jet collection are required for the ECAL bad crystal removal (higgs_dna/tools/EcalBadCalibCrystal_events.py), HiggsDNA crashes without this.
    ],
    "PuppiMET": [
        "pt",
        "phi", # these specific fields from the Jet collection are required for the ECAL bad crystal removal (higgs_dna/tools/EcalBadCalibCrystal_events.py), HiggsDNA crashes without these
    ],
}

SCALARS = [
    "run",
    "luminosityBlock",
    "event",
    # Trigger branches
    "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId",
    "Rho_fixedGridRhoAll", # not sure about its purpose at the moment
]

# MC-only scalar branches
SCALARS_MC = [
    "GenPart"
    "genWeight",
    "Pileup_nTrueInt",
    "Pileup_nPU",
]

WEIGHTS = [
#    "PSWeight",
#    "LHEScaleWeight",
#    "LHEPdfWeight",
]

