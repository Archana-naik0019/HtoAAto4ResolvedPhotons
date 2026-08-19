import awkward as ak
import uproot

from core.reader import NanoReader

reader = NanoReader(
    "root://xrootd-cms.infn.it///store/mc/RunIII2024Summer24NanoAODv15/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v3/2810000/f6d8a890-db4a-4605-99aa-3c30e2a18b19.root"
)

jets = reader.read("Jet")

# Zip all Jet fields into one record
jet = ak.zip({field: jets[field] for field in jets.fields})

with uproot.recreate("test.root") as f:
    tree = f.mktree("Events", {"Jet": jet.type})
    tree.extend({"Jet": jet})

print("Done.")
