import argparse
from h4g_tools.utils.loading import loadPerMass, loadSamples

"""
Only HiggsDNA outputs to check the total number of events. This is only a testing script!
"""


if __name__ == "__main__":
   
    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser")
    parser.add_argument("-y", "--year", type=str, required=True, help="Year to load.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-perMass", "--perMass", type=str, help="Path to files to load per mass.")
    group.add_argument("-input", "--input", type=str, help="Path to files to load all together (combining eras/masses).")

    args = parser.parse_args()
    masses = [x for x in range(15, 65, 5)]

    # Note the custom formatting for eras when handling BDT outputs!
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
    elif args.year == "2023postBPix":
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
    elif args.year == "2024":
        eras={
            "data": ["Run2024C", "Run2024D", "Run2024E", "Run2024F", "Run2024G", "Run2024H", "Run2024I"],
            "bkg": ["RunEvtMix2024C", "RunEvtMix2024D", "RunEvtMix2024E", "RunEvtMix2024F", "RunEvtMix2024G", "RunEvtMix2024H", "RunEvtMix2024I"],
            "signal": [f"Signal_{m}_GeV" for m in masses]
        }
    else:
        assert False, "Must specify year in Scales class!"


    if args.perMass is not None:
        samples = loadPerMass(args.perMass, eras=eras, meta=True, variation="nominal", branches=["pho1_pt"], year=args.year)
        for mass in samples.keys():
            print(samples[mass].schema.metadata[b'initial_events'])

    elif args.input is not None:
        samples = loadSamples(args.input, ["pho1_pt"])
