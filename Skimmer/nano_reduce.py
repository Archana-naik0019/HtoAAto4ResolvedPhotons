import argparse

from core.reducer import NanoReducer
from core.writer import NanoWriter

parser = argparse.ArgumentParser()

parser.add_argument("--input", required=True, help="Input NanoAOD ROOT file")
parser.add_argument("--output", default="skim.root", help="Output skimmed ROOT file")

# New arguments for Lumi Mask and Trigger
parser.add_argument(
    "--is-data",
    action="store_true",
    help="Treat input as Data (enables Lumi Mask and Trigger filtering).",
)
parser.add_argument(
    "--lumimask",
    default=None,
    help="Path to Golden JSON file for lumi mask filtering.",
)
parser.add_argument(
    "--triggers",
    nargs="+",
    default=["HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId"],
    help="List of HLT trigger paths to apply.",
)

# Existing selection flags
parser.add_argument(
    "--no-jet-selection",
    action="store_true",
    help="Disable jet selection.",
)

parser.add_argument(
    "--no-electron-selection",
    action="store_true",
    help="Disable electron selection.",
)

parser.add_argument(
    "--no-muon-selection",
    action="store_true",
    help="Disable muon selection.",
)

parser.add_argument(
    "--no-photon-selection",
    action="store_true",
    help="Disable photon selection.",
)

parser.add_argument(
    "--no-event-selection",
    action="store_true",
    help="Disable event selection.",
)

args = parser.parse_args()

# Configure lumi mask & trigger based on --is-data
apply_lumimask = args.is_data and (args.lumimask is not None)
apply_trigger = args.is_data

store = NanoReducer(
    args.input,
    jet_selection=not args.no_jet_selection,
    electron_selection=not args.no_electron_selection,
    muon_selection=not args.no_muon_selection,
    photon_selection=not args.no_photon_selection,
    event_selection=not args.no_event_selection,
    apply_lumimask=apply_lumimask,
    lumimask_json=args.lumimask,
    apply_trigger=apply_trigger,
    trigger_paths=args.triggers,
).run()

NanoWriter(args.output).write(store)