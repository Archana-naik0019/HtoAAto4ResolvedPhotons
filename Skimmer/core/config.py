MAX_EVENTS_PER_FILE = 500000

COLLECTIONS = [
#    "Jet", # won't need in my analysis
    "Photon",
    "Electron",
#    "Muon", # won't need in my analysis
#    "PuppiMET", # won't need in my analysis
#    "PFMET",
#    "PV",
#    "GenPart",
#    "Flag",
#    "GenVtx",
]

SCALARS = [
    "run",
    "luminosityBlock",
    "event",
#    "Rho_fixedGridRhoFastjetAll",
#    "Rho_fixedGridRhoAll",
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

