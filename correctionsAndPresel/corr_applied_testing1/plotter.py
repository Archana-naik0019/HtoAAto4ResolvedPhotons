#!/usr/bin/env python
"""
plot_corrections.py

Overlay pt and energy distributions (corrected vs. uncorrected)
for the 4 leading photons.
"""

import argparse
import matplotlib
matplotlib.use("Agg")  # no display on lxplus by default
import matplotlib.pyplot as plt
import pandas as pd


def get_parser():
    parser = argparse.ArgumentParser(description="Overlay corrected vs uncorrected photon pt/energy")
    parser.add_argument("--corr", required=True, type=str, help="Path to corrected photons parquet")
    parser.add_argument("--nocorr", required=True, type=str, help="Path to uncorrected photons parquet")
    parser.add_argument("--outdir", default=".", type=str, help="Output directory for plots")
    parser.add_argument("--bins", default=50, type=int, help="Number of histogram bins")
    return parser


def overlay_hist(df_corr, df_nocorr, var, outdir, bins):
    fig, ax = plt.subplots(figsize=(7, 5))

    # Determine common binning range from both datasets combined, for fair comparison
    combined = pd.concat([df_corr[var], df_nocorr[var]])
    lo, hi = combined.min(), combined.max()

    ax.hist(
        df_nocorr[var], bins=bins, range=(lo, hi),
        histtype="step", linewidth=1.8, label="No correction", color="tab:blue",
    )
    ax.hist(
        df_corr[var], bins=bins, range=(lo, hi),
        histtype="step", linewidth=1.8, label="Scale_IJazZ corrected", color="tab:red",
    )

    ax.set_xlabel(var)
    ax.set_ylabel("Events")
    ax.set_title(var)
    ax.legend()
    fig.tight_layout()

    outpath = f"{outdir}/{var}_overlay.png"
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"Saved: {outpath}")


def main():
    args = get_parser().parse_args()

    df_corr = pd.read_parquet(args.corr)
    df_nocorr = pd.read_parquet(args.nocorr)

    # Align on event identifiers so we're comparing the same events
    common_keys = ["run", "luminosityBlock", "event"]
    merged_keys = pd.merge(
        df_corr[common_keys], df_nocorr[common_keys],
        on=common_keys, how="inner",
    )
    print(f"Corrected events: {len(df_corr)}, Uncorrected events: {len(df_nocorr)}, "
          f"Common events: {len(merged_keys)}")

    for i in range(1, 5):
        for quantity in ["pt", "energy"]:
            var = f"pho{i}_{quantity}"
            if var not in df_corr.columns or var not in df_nocorr.columns:
                print(f"Skipping {var}: not found in one of the dataframes")
                continue
            overlay_hist(df_corr, df_nocorr, var, args.outdir, args.bins)


if __name__ == "__main__":
    main()
