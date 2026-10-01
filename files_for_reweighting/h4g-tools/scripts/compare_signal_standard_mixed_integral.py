import os
import argparse
from h4g_tools.utils.scales import Scales as ScalesCls
from h4g_tools.utils.loading import loadPerMass, fillHist
import gc


if __name__ == "__main__":

    gc.disable()

    parser = argparse.ArgumentParser()
    parser.add_argument("-std", "--standard", type=str, required=True, help="Standard signal samples.")
    parser.add_argument("-mix", "--mixed", type=str, required=True, help="Mixed signal samples.")
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    args = parser.parse_args()
    
    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Scales instance setup
    Scales = ScalesCls(args.year) 
    masses = [x for x in range(15, 65, 5)]
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
            "data": ["Run2023D"],
            "bkg": ["RunEvtMix2023D"],
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

    # Load samples
    print("Loading signal samples...")
    standard = loadPerMass(args.standard, eras = eras, year = args.year, branches=["BDT_score", "weight"])
    mixed = loadPerMass(args.mixed, eras = eras, year = args.year, branches=["BDT_score", "weight"])

    sorted_samples = sorted(standard.keys())
    for idx, (mass, score) in enumerate(zip(sorted_samples, [0.9500, 0.9550, 0.9650, 0.9800, 0.9800, 0.9800, 0.9800, 0.9800, 0.9800, 0.9800])):
        std = standard[mass]
        mix = mixed[mass]

        std_mass = fillHist("BDT_score", std, [0.0, 1.0], binsScale=200, name=f"std_BDT_{mass[7:13]}", normalize=False)
        std_mass.Scale(Scales.lumi_fb * 1.0 / sum(std.weight))
        mix_mass = fillHist("BDT_score", mix, [0.0, 1.0], binsScale=200, name=f"mix_BDT_{mass[7:13]}", normalize=False)
        mix_mass.Scale(Scales.lumi_fb * 1.0 / sum(mix.weight))

        std_int = std_mass.Integral(std_mass.FindBin(score), std_mass.GetNbinsX())
        mix_int = mix_mass.Integral(mix_mass.FindBin(score), mix_mass.GetNbinsX())
        ratio_std_mix = std_int / mix_int
        ratio_mix_std = mix_int / std_int

        print(mass, score, f"standard: {std_int}", f"mixed: {mix_int}", f"ratio std/mix: {ratio_std_mix}", f"ratio mix/std: {ratio_mix_std}")
        print()


