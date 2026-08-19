import uproot

f = uproot.open("root://xrootd-cms.infn.it///store/mc/RunIII2024Summer24NanoAODv15/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v3/2810000/f6d8a890-db4a-4605-99aa-3c30e2a18b19.root")

tree = f["Events"]

for b in tree.keys():
    if "Rho" in b or "rho" in b:
        print(b)
