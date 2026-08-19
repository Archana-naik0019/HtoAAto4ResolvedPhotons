import argparse

from core.reducer import NanoReducer
from core.writer import NanoWriter

parser = argparse.ArgumentParser()

parser.add_argument("--input", required=True)
parser.add_argument("--output", default="skim.root")

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

store = NanoReducer(
    args.input,
    jet_selection=not args.no_jet_selection,
    electron_selection=not args.no_electron_selection,
    muon_selection=not args.no_muon_selection,
    photon_selection=not args.no_photon_selection,
    event_selection=not args.no_event_selection,
).run()

NanoWriter(args.output).write(store)