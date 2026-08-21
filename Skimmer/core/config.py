MAX_EVENTS_PER_FILE = 500000

COLLECTIONS = [
    "Jet", # won't need in my analysis
#    "Jet_pt",
#    "Jet_eta",
#    "Jet_phi",
#    "Jet_neEmEF",
#    "Jet_chEmEF",
    "Photon",
    "Electron",
#    "Muon", # won't need in my analysis
    "PuppiMET", # won't need in my analysis
#    "PuppiMET_pt",
#    "PuppiMET_phi",
#    "PFMET",
    "PV",
#    "GenPart",
#    "Flag",
#    "GenVtx",
]

KEEP_FIELDS = {
    "Jet": [
        "pt",
        "eta",
        "phi",
        "neEmEF",
        "chEmEF",
    ],
    "PuppiMET": [
        "pt",
        "phi",
    ],
}

SCALARS = [
    "run",
    "luminosityBlock",
    "event",
    # Trigger branches
    "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId",
    "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55",
#    "Rho_fixedGridRhoFastjetAll",
    "Rho_fixedGridRhoAll",
]

# MC-only scalar branches
SCALARS_MC = [
    "genWeight",
    "Pileup_nTrueInt",
    "Pileup_nPU",
]

WEIGHTS = [
#    "PSWeight",
#    "LHEScaleWeight",
#    "LHEPdfWeight",
]

