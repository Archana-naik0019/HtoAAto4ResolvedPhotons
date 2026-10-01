from h4g_tools.utils.loading import loadPerMass, fillHist
from h4g_tools.bdt.event_selection_bdt import reduce_dataframe, predict_BDT
from h4g_tools.utils.scales import Scales as ScalesCls
from pathlib import Path
import xgboost
import os
import time
from copy import deepcopy

if __name__ == "__main__":

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Scales instance setup
    year = "2024"
    masses = [x for x in range(15, 65, 5)]
    eras={
        "data": ["Run2024C", "Run2024D", "Run2024E", "Run2024F", "Run2024G", "Run2024H", "Run2024I"],
        "bkg": ["RunEvtMix2024C", "RunEvtMix2024D", "RunEvtMix2024E", "RunEvtMix2024F", "RunEvtMix2024G", "RunEvtMix2024H", "RunEvtMix2024I"],
        "signal": [f"Signal_{m}_GeV" for m in masses]
    }

    # Handle path setup after checking year
    Scales = ScalesCls(year) 

    # Load samples (nominal only by default)
    print("Loading signal samples...")
    sig_input_path = "signal-input-path-to-be-varied"
    samples = loadPerMass(sig_input_path)
    dump = "outputs_signal2024_photonMVAID_variations/"

    # 8 variations since we're doing it one-by-one
    # pho1 +/- variation%
    # pho2 +/- variation%
    # pho3 +/- variation%
    # pho4 +/- variation%

    # Variation multipliers
    #variation = 0.1
    #variation = 0.2
    variation = 0.3
    var_up   = 1.0 + variation
    var_down = 1.0 - variation

    # Load BDT model
    bdt = xgboost.XGBClassifier()
    bdt.load_model(f"{cwd}/../models/ES_BDT_2024_updatedPresel.xgb")

    # Do one variation for one photon at a time per full sample
    current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    vrs = [var_up, var_down]
    for idx, var in enumerate(vrs):
        print(f"\nProcessing {'up' if var == var_up else 'down'} variation:")
        for pho in range(1,5):
            print(f"~ Photon {pho}")
            for mass in sorted(samples.keys()):
                print(f"~~~ Mass: {mass.replace('Signal_','').replace('_',' ')}")
                s = deepcopy(samples[mass])

                # Replace existing variable with variation so it can be used by BDT
                gt = s[f"pho{pho}_mvaID"] >= 0.0
                lt = s[f"pho{pho}_mvaID"] < 0.0
                s[f"pho{pho}_mvaID"][gt] = s[f"pho{pho}_mvaID"][gt] * var
                s[f"pho{pho}_mvaID"][lt] = s[f"pho{pho}_mvaID"][lt] * vrs[(idx+1)%2]
                # > 0.0 -> can do +/- 10%
                # < 0.0 -> need to do -0.1 + 10% = -0.9 = -0.1 * (0.9 NOT 1.1)

                # Remove MVA scores < -1 and > 1
                s[f"pho{pho}_mvaID"].loc[s[f"pho{pho}_mvaID"] > 1.0] = 1.0
                s[f"pho{pho}_mvaID"].loc[s[f"pho{pho}_mvaID"] < -1.0] = -1.0

                # Don't overwrite existing DataFrame -> just dump current

                # Add BDT evaluation as new column
                s["BDT_score"] = predict_BDT(bdt, reduce_dataframe(s))

                # Will need to dump varied samples into a new folder (maybe just variation, maybe an entirely new directory?)
                out = f"pho{pho}_mvaID_varUp" if var == var_up else f"pho{pho}_mvaID_varDown"
                outpath = os.path.join(*[dump, out, mass, "mvaID_Up" if var == var_up else "mvaID_Down"])
                filepath = f"merged_{current_time}.parquet"

                Path(outpath).mkdir(parents=True, exist_ok=True)
                s.to_parquet(os.path.join(*[outpath, filepath]))


