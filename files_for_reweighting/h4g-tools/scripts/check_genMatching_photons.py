import os
import argparse
import json
from ROOT import TLorentzVector, gErrorIgnoreLevel, kFatal
import numpy as np
from h4g_tools.utils.scales import Scales as ScalesCls
from h4g_tools.utils.loading import pathSetup, loadPerMass, fillHist

gErrorIgnoreLevel = kFatal

def applyCuts(
    samples: dict,
    valid: list
) -> dict:

    pass


def fillVectors(
    vecs: list,
    inputs: list
) -> bool:
    
    assert len(vecs) == len(inputs), (len(vecs), len(inputs))
    any_nan = False

    for jdx, (inp, vec) in enumerate(zip(inputs, vecs)):
        vec.SetPtEtaPhiM(
            inp[0],
            inp[1],
            inp[2],
            inp[3],
        )

        # All four fields would be nan simultaneously by construction, thus save on compute by only checking one.
        if np.isnan(inp[0]):
            any_nan = True

    return any_nan


def checkGenMatching(
    samples: dict,
    dR_min: float,
    year: str,
    saveJSON: bool
) -> None:

    valid = []

    # Do setup for vectors, then pass them to function for filling
    # 4vectors for all 4 generator photons
    gen_pho1 = TLorentzVector(0,0,0,0)
    gen_pho2 = TLorentzVector(0,0,0,0)
    gen_pho3 = TLorentzVector(0,0,0,0)
    gen_pho4 = TLorentzVector(0,0,0,0)

    # 4vectors for all 4 RECO photons
    pho1 = TLorentzVector(0,0,0,0)
    pho2 = TLorentzVector(0,0,0,0)
    pho3 = TLorentzVector(0,0,0,0)
    pho4 = TLorentzVector(0,0,0,0)

    # Do gen-matching checks to see success rate per mass. Can combine it all later if desired.
    eff = {"dR_min": dR_min}
    samples_keys = sorted(samples.keys()) 
    for k in samples_keys:
        s = samples[k]
        mass = k.replace("Signal_","")
        mass = mass[:mass.find("GeV")+3]
        print(f"Processing {mass.replace('_',' ')}...")

        # Match both = dR match for both pseudoscalars
        # No both match = dR is too big for a match for both pseudoscalars
        # No a1 match = dR is too big for a match of a1
        # No a2 match = dR is too big for a match of a2
        # No gen Idx = not possible to compare with RECO since there is no corresponding gen-level object
        eff.update({mass: {"total_events": 0, "match_both": 0, "match_a1": 0, "match_a2": 0, "no_both_match": 0, "no_a1": 0, "no_a2": 0, "no_genIdx": 0}})

        # Load gen_pho1
        gen_pho1_pt = s.gen_pho1_pt.to_numpy()
        gen_pho1_eta = s.gen_pho1_eta.to_numpy()
        gen_pho1_phi = s.gen_pho1_phi.to_numpy()
        gen_pho1_mass = s.gen_pho1_mass.to_numpy()

        # Load gen_pho2
        gen_pho2_pt = s.gen_pho2_pt.to_numpy()
        gen_pho2_eta = s.gen_pho2_eta.to_numpy()
        gen_pho2_phi = s.gen_pho2_phi.to_numpy()
        gen_pho2_mass = s.gen_pho2_mass.to_numpy()

        # Load gen_pho3
        gen_pho3_pt = s.gen_pho3_pt.to_numpy()
        gen_pho3_eta = s.gen_pho3_eta.to_numpy()
        gen_pho3_phi = s.gen_pho3_phi.to_numpy()
        gen_pho3_mass = s.gen_pho3_mass.to_numpy()

        # Load gen_pho4
        gen_pho4_pt = s.gen_pho4_pt.to_numpy()
        gen_pho4_eta = s.gen_pho4_eta.to_numpy()
        gen_pho4_phi = s.gen_pho4_phi.to_numpy()
        gen_pho4_mass = s.gen_pho4_mass.to_numpy()

        if not applyCuts:
            pho1_str = "LeadPs_leading_pho_"
            pho2_str = "LeadPs_subleading_pho_"
            pho3_str = "SubleadPs_leading_pho_"
            pho4_str = "SubleadPs_subleading_pho_"
        else:
            pho1_str = "pho1_"
            pho2_str = "pho2_"
            pho3_str = "pho3_"
            pho4_str = "pho4_"

        # Load RECO pho1
        pho1_pt = s[f"{pho1_str}pt"].to_numpy()
        pho1_eta = s[f"{pho1_str}eta"].to_numpy()
        pho1_phi = s[f"{pho1_str}phi"].to_numpy()
        pho1_mass = s[f"{pho1_str}mass"].to_numpy()

        # Load RECO pho2
        pho2_pt = s[f"{pho2_str}pt"].to_numpy()
        pho2_eta = s[f"{pho2_str}eta"].to_numpy()
        pho2_phi = s[f"{pho2_str}phi"].to_numpy()
        pho2_mass = s[f"{pho2_str}mass"].to_numpy()

        # Load RECO pho3
        pho3_pt = s[f"{pho3_str}pt"].to_numpy()
        pho3_eta = s[f"{pho3_str}eta"].to_numpy()
        pho3_phi = s[f"{pho3_str}phi"].to_numpy()
        pho3_mass = s[f"{pho3_str}mass"].to_numpy()

        # Load RECO pho4
        pho4_pt = s[f"{pho4_str}pt"].to_numpy()
        pho4_eta = s[f"{pho4_str}eta"].to_numpy()
        pho4_phi = s[f"{pho4_str}phi"].to_numpy()
        pho4_mass = s[f"{pho4_str}mass"].to_numpy()

        for idx in range(len(s)):
            #print(f"~~~~ Event {idx} ~~~~")
            nan = False

            # Fill ROOT vectors for RECO photons
            nan = fillVectors(
                [pho1, pho2, pho3, pho4],
                [
                    [pho1_pt[idx], pho1_eta[idx], pho1_phi[idx], pho1_mass[idx]],
                    [pho2_pt[idx], pho2_eta[idx], pho2_phi[idx], pho2_mass[idx]],
                    [pho3_pt[idx], pho3_eta[idx], pho3_phi[idx], pho3_mass[idx]],
                    [pho4_pt[idx], pho4_eta[idx], pho4_phi[idx], pho4_mass[idx]],
                ]
            )

            # Fill ROOT vectors for generator photons
            nan = nan or fillVectors(
                [gen_pho1, gen_pho2, gen_pho3, gen_pho4],
                [
                    [gen_pho1_pt[idx], gen_pho1_eta[idx], gen_pho1_phi[idx], gen_pho1_mass[idx]],
                    [gen_pho2_pt[idx], gen_pho2_eta[idx], gen_pho2_phi[idx], gen_pho2_mass[idx]],
                    [gen_pho3_pt[idx], gen_pho3_eta[idx], gen_pho3_phi[idx], gen_pho3_mass[idx]],
                    [gen_pho4_pt[idx], gen_pho4_eta[idx], gen_pho4_phi[idx], gen_pho4_mass[idx]],
                ]
            )

            # Only look for valid events if not saving JSON (intentionally)
            if not saveJSON:
                print(gen_pho1.Pt(), pho1.Pt(), gen_pho2.Pt(), pho2.Pt(), gen_pho3.Pt(), pho3.Pt(), gen_pho4.Pt(), pho4.Pt())
                continue
                if gen_pho1.DeltaR(pho1) < 0.1 and gen_pho2.DeltaR(pho2) < 0.1 and gen_pho3.DeltaR(pho3) < 0.1 and gen_pho4.DeltaR(pho4) < 0.1:
                    valid.append(idx)
                    continue

            # Check dR per event between generator and RECO photons
            # match a generator and RECO photon
            # Check if pseudoscalars are built correctly by combining the matching gen-level photons into pseudoscalar objects knowing which Ps the RECO objects came from
            # i.e., pho1 + pho2 = ps1 and pho3 + pho4 = ps2
            # the total success rate will be the reconstruction efficiency in signal MC. works per mass point or all together

            gen_ps11 = gen_pho1 + gen_pho2
            gen_ps12 = gen_pho3 + gen_pho4
            gen_ps21 = gen_pho1 + gen_pho3
            gen_ps22 = gen_pho2 + gen_pho4
            gen_ps31 = gen_pho1 + gen_pho4
            gen_ps32 = gen_pho2 + gen_pho4

            # Not sure which photons make which pseudoscalars, so need to do all permutations
            if gen_ps11.Pt() >= gen_ps12.Pt():
                gen1_a1 = gen_ps11
                gen1_a2 = gen_ps12
            else:
                gen1_a1 = gen_ps12
                gen1_a2 = gen_ps11
            if gen_ps21.Pt() >= gen_ps22.Pt():
                gen2_a1 = gen_ps21
                gen2_a2 = gen_ps22
            else:
                gen2_a1 = gen_ps22
                gen2_a2 = gen_ps21
            if gen_ps31.Pt() >= gen_ps32.Pt():
                gen3_a1 = gen_ps31
                gen3_a2 = gen_ps32
            else:
                gen3_a1 = gen_ps32
                gen3_a2 = gen_ps31
            
            # RECO ordering is set by construction
            reco_ps1 = pho1 + pho2
            reco_ps2 = pho3 + pho4

            # NaNs come from particles which have no corresponding/matching generator object (e.g. photons from Events that cannot be matched to the a->gg decay)
            # See under "IMPORTANT" about "stable": https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookGenParticleCandidate#GenPCand
            # and under GenPart_status here: https://cms-xpog.docs.cern.ch/autoDoc/NanoAODv15/2024/doc_TTH-Hto2G_Par-M-125_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2.html

            a1_check = reco_ps1.DeltaR(gen1_a1) < dR_min or reco_ps1.DeltaR(gen2_a1) < dR_min or reco_ps1.DeltaR(gen3_a1) < dR_min
            a2_check = reco_ps2.DeltaR(gen1_a2) < dR_min or reco_ps2.DeltaR(gen2_a2) < dR_min or reco_ps2.DeltaR(gen3_a2) < dR_min

            eff[mass]["total_events"] += 1

            # Check NaNs
            if nan:
                eff[mass]["no_genIdx"] += 1
                # Don't double count NaN and no match!
                continue

            # Check a1
            if a1_check:
                eff[mass]["match_a1"] += 1
            else:
                eff[mass]["no_a1"] += 1

            # Check a2
            if a2_check:
                eff[mass]["match_a2"] += 1
            else:
                eff[mass]["no_a2"] += 1

            # Check both a1 and a2 together
            if a1_check and a2_check:
                eff[mass]["match_both"] += 1
            else:
                eff[mass]["no_both_match"] += 1

    # Save reconstruction efficiency to JSON
    if not saveJSON:
        with open(f"{cwd}/gen_matching_check_{year}_dR{str(dR_min).replace('.','p')}.json", "w") as f:
            json.dump(eff, f, indent=4)
    else:
        return valid


if __name__ == "__main__":
    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser", conflict_handler="resolve")

    # Input arguments
    parser.add_argument("-i", "--inputs", help="Directory to import signal samples from.", type=str, default=None, required=True)
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    parser.add_argument("-dR", "--dR_min", type=float, required=True, help="dR minimum for gen/reco matching.")
    parser.add_argument("-cuts", "--applyCuts", action="store_true", help="Apply h4g selections to samples after matching.")
    args = parser.parse_args()

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Scales instance setup
    Scales = ScalesCls(args.year) 
    masses = [x for x in range(15, 65, 5)]
    if args.year == "2022preEE":
        eras={
            "data": ["Run2022C","Run2022D"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D"],
            "signal": [f"Signal_{m}_GeV_preEE" for m in masses]
        }
    elif args.year == "2022postEE":
        eras={
            "data": ["Run2022E","Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_postEE" for m in masses]
        }
    elif args.year == "2022":
        eras={
            "data": ["Run2022C","Run2022D", "Run2022E","Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D", "RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_{prepost}" for m in masses for prepost in ["preEE", "postEE"]]
        }
    elif args.year == "2024":
        eras={
            "data": ["Run2024C", "Run2024D", "Run2024E", "Run2024F", "Run2024G", "Run2024H", "Run2024I"],
            "bkg": ["RunEvtMix2024C", "RunEvtMix2024D", "RunEvtMix2024E", "RunEvtMix2024F", "RunEvtMix2024G", "RunEvtMix2024H", "RunEvtMix2024I"],
            "signal": [f"Signal_{m}_GeV" for m in masses]
        }
    else:
        assert False, "Must specify year in Scales class!"

    # Handle path setup after checking year
    pathSetup(args.year)

    # Load signal MC samples per mass point
    samples = loadPerMass(args.inputs, eras = eras, year = args.year, era_type_arg = "signal", skipInterMass=True)

    valid = checkGenMatching(samples, args.dR_min, args.year, not args.applyCuts)
    print(valid)
    if args.applyCuts:
        # Do h4g cuts and gen matching ensuring that the reco photons having a matching gen photon first
        samples = applyCuts(samples, valid)
        checkGenMatching(samples, args.dR_min, args.year, args.applyCuts)
