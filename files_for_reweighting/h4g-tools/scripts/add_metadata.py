import argparse
import os
import pyarrow
import json
from pathlib import Path
from h4g_tools.utils.loading import loadPerMass, loadSamples


"""
Example:
python3 add_metadata.py -meta /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2022_SaS_Pileup_NFs_Material_FNUF_eVetoSF_PreselSF_TriggerSF_08Dec2025/ -y 2022 -xgb ES_BDT_2022 -bdt /users/scastel2/h4g-tools/scripts/cats/2022/ES_BDT_2022/bdt_2022_cats1.txt

Can run signal and background independently!

NOTE: Use the output of apply_categories.py for the -bdt flag!
"""



if __name__ == "__main__":
   
    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser")
    parser.add_argument("-meta", "--metadata", type=str, required=False, help="Path to files with metadata included.")
    parser.add_argument("-bkg_meta", "--bkg_metadata", type=str, required=False, help="Path to files with metadata included for background.")
    parser.add_argument("-y", "--year", type=str, required=True, help="Year to load.")
    parser.add_argument("-xgb", "--xgb_name", type=str, required=True, help="XGBoost model name.")
    parser.add_argument("-bdt", "--bdt_path", required=True, type=str, help="Path to file containing efficiency of BDT cuts. Assumes 1 category.")
    parser.add_argument("-t", "--test", action="store_true", required=False, help="Only run nominal to test metadata outputs.")
    args = parser.parse_args()

    masses = [x for x in range(15, 65, 5)]

    # Note the custom formatting for eras when handling BDT outputs!
    if args.year == "2022preEE":
        eras={
            "data": ["Run2022C", "Run2022D"],
            "bkg": [f"{m}_GeV" for m in masses],
            "signal": [f"2022preEE/{args.xgb_name}/signal/{m}_GeV" for m in masses]
        }
        eras_meta={
            "data": ["Run2022C", "Run2022D"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D"],
            "signal": [f"Signal_{m}_GeV_preEE" for m in masses]
        }
    elif args.year == "2022postEE":
        eras={
            "data": ["Run2022E", "Run2022F", "Run2022G"],
            "bkg": [f"{m}_GeV" for m in masses],
            "signal": [f"2022postEE/{args.xgb_name}/signal/{m}_GeV" for m in masses]
        }
        eras_meta={
            "data": ["Run2022E", "Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_postEE" for m in masses]
        }
    elif args.year == "2022":
        eras={
            "data": ["Run2022C", "Run2022D", "Run2022E", "Run2022F", "Run2022G"],
            "bkg": [f"{m}_GeV" for m in masses],
            "signal": [f"{args.year}/{args.xgb_name}/signal/{m}_GeV" for m in masses]
        }
        eras_meta={
            "data": ["Run2022C", "Run2022D", "Run2022E", "Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D", "RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_{prepost}" for m in masses for prepost in ["preEE", "postEE"]]
        }
    elif args.year == "2023preBPix":
        eras={
            "data": ["Run2023C"],
            "bkg": [f"{m}_GeV" for m in masses],
            "signal": [f"2023preBPix/{args.xgb_name}/signal/{m}_GeV" for m in masses]
        }
        eras_meta={
            "data": ["Run2023C"],
            "bkg": ["RunEvtMix2023C"],
            "signal": [f"Signal_{m}_GeV_preBPix" for m in masses]
        }
    elif args.year == "2023postBPix":
        eras={
            "data": ["Run2023D"],
            "bkg": [f"{m}_GeV" for m in masses],
            "signal": [f"2023postBPix/{args.xgb_name}/signal/{m}_GeV" for m in masses]
        }
        eras_meta={
            "data": ["Run2022D"],
            "bkg": ["RunEvtMix2022D"],
            "signal": [f"Signal_{m}_GeV_postBPix" for m in masses]
        }
    elif args.year == "2023":
        eras={
            "data": ["Run2023C", "Run2023D"],
            "bkg": [f"{m}_GeV" for m in masses],
            "signal": [f"2023/{args.xgb_name}/signal/{m}_GeV" for m in masses]
        }
        eras_meta={
            "data": ["Run2023C", "Run2023D"],
            "bkg": ["RunEvtMix2023C", "RunEvtMix2023D"],
            "signal": [f"Signal_{m}_GeV_{prepost}" for m in masses for prepost in ["preBPix", "postBPix"]]
        }
    elif args.year == "2024":
        eras={
            "data": ["Run2024C", "Run2024D", "Run2024E", "Run2024F", "Run2024G", "Run2024H", "Run2024I"],
            "bkg": [f"{m}_GeV" for m in masses],
            "signal": [f"{args.year}/{args.xgb_name}/signal/{m}_GeV" for m in masses]
        }
        eras_meta={
            "data": ["Run2024C", "Run2024D", "Run2024E", "Run2024F", "Run2024G", "Run2024H", "Run2024I"],
            "bkg": ["RunEvtMix2024C", "RunEvtMix2024D", "RunEvtMix2024E", "RunEvtMix2024F", "RunEvtMix2024G", "RunEvtMix2024H", "RunEvtMix2024I"],
            "signal": [f"Signal_{m}_GeV" for m in masses]
        }
    else:
        assert False, "Must specify year in Scales class!"

    print()

    metadata_path = args.metadata
    bkg_metadata_path = args.bkg_metadata
    sig_path = f"/cms/cephfs/data/store/user/castells/outputs/BDT/{args.year}/{args.xgb_name}/signal/"
    bkg_path = f"/cms/cephfs/data/store/user/castells/outputs/BDT/{args.year}/{args.xgb_name}/background/"

    # Dictionary setup
    samples = {}
    bkg_samples = {}
    metadata = {}
    bkg_metadata = {}
    metadata_signal_export = {}

    # Load each variation for signal
    if not args.test and args.metadata is not None:
        variations = os.listdir(os.path.join(*(sig_path, os.listdir(sig_path)[0])))
    elif args.test:
        variations = ["nominal"]
    elif args.metadata is None:
        variations = ["nominal"]
    print(variations)
    for variation in variations:
        if args.metadata is not None:
            print(f"=== Loading variation: {variation} ===")
            print("~~~ Loading signal MC from BDT outputs ~~~")
            samples.update({variation: {}})
            samples[variation].update(loadPerMass(sig_path, eras=eras, variation=variation, bdt_outputs=True, year=args.year))

            print("\n~~~ Loading signal MC metadata from HiggsDNA outputs ~~~")
            metadata.update({variation: {}})
            metadata[variation].update(loadPerMass(metadata_path, eras=eras_meta, variation=variation, meta=True, branches=["pho1_pt"], year=args.year))

        # Load bkg samples
        if variation == "nominal" and args.bkg_metadata is not None:
            print("\n~~~ Loading event mixed background from BDT outputs ~~~")
            bkg_samples.update({variation: {}})
            bkg_samples[variation].update(loadPerMass(bkg_path, eras=eras, variation=variation, bdt_outputs=True, year=args.year))

            print("\n~~~ Loading event mixed background metadata from HiggsDNA outputs ~~~")
            bkg_metadata.update({variation: loadSamples(bkg_metadata_path, eras=eras_meta, variation=variation, meta=True, branches=["pho1_pt"])})

        assert samples.keys() == metadata.keys()
        if args.metadata is not None:
            print("\n")

    with open(args.bdt_path, "r") as f:
        bdt_eff = json.load(f)

    # Dump BDT outputs with original metadata included
    for variation in variations:
        print(f"\n~~~ Processing variation: {variation} ~~~")
        metadata_signal_export.update({variation: {}})
        if args.metadata is not None:
            for mass in [str(m) for m in range(15,65,5)]:
                mass_new = mass + "_GeV"
                outpath = os.path.join(*[sig_path, mass_new, variation])
                print(f"Running over mass {mass_new.replace('_',' ')}...", end="\t")

                # BDT cut checks
                bdt_cut = bdt_eff[mass_new.replace("Signal_","")]["cut0"][0]
                sum_of_weights = sum(samples[variation][mass_new][samples[variation][mass_new].BDT_score > bdt_cut]["weight_central"])

                # Add metadata to samples with BDT included
                new_meta = metadata[variation][mass].schema.metadata
                new_meta.update({b'bdt': str(len(samples[variation][mass_new][samples[variation][mass_new].BDT_score > bdt_cut])).encode("utf-8")})
                #new_meta[b'sum_weight_central'] = str(sum_of_weights).encode("utf-8")  # Update this to the sum of weights AFTER bdt cuts applied
                samples[variation][mass_new] = pyarrow.Table.from_pandas(samples[variation][mass_new]).replace_schema_metadata(new_meta)

                str_new_meta = {}
                for k in new_meta:
                    str_new_meta.update({str(k.decode("utf-8")): str(new_meta[k].decode("utf-8"))})
                metadata_signal_export[variation].update({mass: str_new_meta})

                # Corresponding variation path
                variation_path = os.path.join(*[sig_path, mass_new, variation])

                assert len(os.listdir(variation_path)) == 1, os.listdir(variation_path)
                file_path = os.path.join(*[variation_path, os.listdir(variation_path)[0]])

                # Write pyarrow table to parquet file
                if not args.test:
                    pyarrow.parquet.write_table(samples[variation][mass_new], file_path)
                    print(f"Wrote {variation} for {mass_new} to file.")

        if args.bkg_metadata is not None and variation == "nominal":
            for mass in [str(m) for m in range(15,65,5)]:
                mass_new = mass + "_GeV"
                outpath = os.path.join(*[bkg_path, mass_new, variation])

                # BDT cut checks
                bdt_cut = bdt_eff[mass_new.replace("Signal_","")]["cut0"][0]
                sum_of_weights = float(len(bkg_samples[variation][mass_new][bkg_samples[variation][mass_new].BDT_score > bdt_cut]))

                # Add metadata to samples with BDT included
                new_meta = bkg_metadata[variation].schema.metadata
                new_meta.update({b'bdt': str(sum_of_weights).encode("utf-8")})
                # sum_weight_central does not exist for bkg, so no need to handle it here
                bkg_samples[variation][mass_new] = pyarrow.Table.from_pandas(bkg_samples[variation][mass_new]).replace_schema_metadata(new_meta)

                # Corresponding variation path
                variation_path = os.path.join(*[bkg_path, mass_new, variation])

                assert len(os.listdir(variation_path)) == 1, os.listdir(variation_path)
                file_path = os.path.join(*[variation_path, os.listdir(variation_path)[0]])

                # Write pyarrow table to parquet file
                if not args.test:
                    spaces = len(f"Running over mass {mass_new.replace('_',' ')}...")
                    pyarrow.parquet.write_table(bkg_samples[variation][mass_new], file_path)
                    print((" " * spaces) + f"Wrote {variation} for {mass_new} background to file.")


    with open(f"meta_signalMC_{args.year}.json", "w") as f:
        json.dump(metadata_signal_export, f, indent=4)

    print("\n")
