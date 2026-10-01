import argparse
import os
import json
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

"""
python3 compare_ams_vs_bdt.py -i1 <path-to-json-for-year1> -i2 <path-to-json-for-year2> -y1 <year1> -y2 <year2>

"""


# Get current directory
cwd = os.path.dirname(os.path.abspath(__file__))


def compareAMS(
    args: argparse.ArgumentParser,
) -> None:

    # Path setup
    if not os.path.exists(os.path.join(cwd, "superimposed_2022_2024_ams_vs_bdt")):
        Path(os.path.join(cwd, "superimposed_2022_2024_ams_vs_bdt")).mkdir(parents=True, exist_ok=True)

    # Define paths for comparison plots
    masses = list(range(15,65,5))
    paths = [os.path.join(cwd, "superimposed_2022_2024_ams_vs_bdt", f"ams_vs_bdt_{args.year1}_{args.year2}_{mass}.png") for mass in masses]

    # Load AMS vs BDT values for each year
    y1 = {}
    y2 = {}
    with open(args.input1, "r") as f:
        y1 = json.load(f)
    with open(args.input2, "r") as f:
        y2 = json.load(f)
    assert y1.keys() == y2.keys()

    # Make plots per mass point
    for mass in y1.keys():
        # Variables for y1
        cut1 = y1[mass]["cut"]
        ams1 = y1[mass]["ams"]
        invalid_cut1 = y1[mass]["invalid_cut"]
        invalid1 = y1[mass]["invalid"]
        max_ams1 = y1[mass]["max_ams"]
        max_ams_noRes1 = y1[mass]["max_ams_noRestrictions"]
        max_cut1 = y1[mass]["axline"]
        data_sb1 = y1[mass]["data_sb"]
        signal_sr1= y1[mass]["signal_sr"]
        bkg_sr1 = y1[mass]["bkg_sr"]

        # Variables for y2
        cut2 = y2[mass]["cut"]
        ams2 = y2[mass]["ams"]
        invalid_cut2 = y2[mass]["invalid_cut"]
        invalid2 = y2[mass]["invalid"]
        max_ams2 = y2[mass]["max_ams"]
        max_ams_noRes2 = y2[mass]["max_ams_noRestrictions"]
        max_cut2 = y2[mass]["axline"]
        data_sb2 = y2[mass]["data_sb"]
        signal_sr2= y2[mass]["signal_sr"]
        bkg_sr2 = y2[mass]["bkg_sr"]

        # Print max AMS and max AMS without restrictions
        """
        print(f"Mass: {mass.replace('_', ' ')}")
        print(f"~~ {args.year1} ~~")
        print(f"   Max AMS:       {max_ams1:.2f}")
        print(f"   Max AMS (all): {max_ams_noRes1:.2f}")
        print(f"~~ {args.year2} ~~")
        print(f"   Max AMS:       {max_ams2:.2f}")
        print(f"   Max AMS (all): {max_ams_noRes2:.2f}")
        """

        # Plot AMS vs BDT cut
        fig, ax1 = plt.subplots(figsize = (8,8))
        title = f"AMS vs BDT Cut for m_{{a}} = {mass.replace('_',' ')}"
        ax1.scatter(cut1, ams1, s=5, c="b", label=args.year1)
        ax1.scatter(invalid_cut1, invalid1, s=5, c="r", label=args.year1)
        ax1.axvline(max_cut1, color="g", ls=":", label=f"{args.year1} Optimal Cut")
        ax1.scatter(cut2, ams2, s=5, c="tab:cyan", marker="s", label=args.year2)
        ax1.scatter(invalid_cut2, invalid2, s=5, c="tab:orange", marker="s", label=args.year2)
        ax1.axvline(max_cut2, color="k", ls="--", label=f"{args.year2} Optimal Cut")
        ax1.set_xlabel("Cut on BDT Score", fontsize=14)
        ax1.set_ylabel("AMS", fontsize=14)
        ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
        ax1.set_ylim(top=39.0, bottom=2.0)
        ax1.legend(fontsize=14)
        
        # Add text to plot
        ax1.text(0.165, 28.0, "BDT Cut", weight="bold", fontsize=13)
        ax1.text(0.135, 27.0, f"{args.year1}: {max_cut1:.4f}", fontsize=13)
        ax1.text(0.135, 26.0, f"{args.year2}: {max_cut2:.4f}", fontsize=13)
        ax1.text(0.13, 24.0, "Optimal AMS", weight="bold", fontsize=13)
        ax1.text(0.15, 23.0, f"{args.year1}: {max_ams1:.2f}", fontsize=13)
        ax1.text(0.15, 22.0, f"{args.year2}: {max_ams2:.2f}", fontsize=13)
        ax1.text(0.10, 20.0, "Max AMS (no data min)", weight="bold", fontsize=13)
        ax1.text(0.15, 19.0, f"{args.year1}: {max_ams_noRes1:.2f}", fontsize=13)
        ax1.text(0.15, 18.0, f"{args.year2}: {max_ams_noRes2:.2f}", fontsize=13)
        ax1.text(0.075, 16.0, "Data in SB at Optimal AMS", weight="bold", fontsize=13)
        ax1.text(0.15, 15.0, f"{args.year1}: {data_sb1:.2f}", fontsize=13)
        ax1.text(0.15, 14.0, f"{args.year2}: {data_sb2:.2f}", fontsize=13)

        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.tight_layout()
        if not os.path.exists(f"{cwd}/plots/BDT/ams_vs_bdt_{args.year1}_{args.year2}/"):
            Path(f"{cwd}/plots/BDT/ams_vs_bdt_{args.year1}_{args.year2}/").mkdir(parents=True, exist_ok=True)
        plt.savefig(f"{cwd}/plots/BDT/ams_vs_bdt_{args.year1}_{args.year2}/ams_vs_bdt_{args.year1}_{args.year2}_{mass}.png")
        print(f"Saved plot ams_vs_bdt_{args.year1}_{args.year2}/ams_vs_bdt_{args.year1}_{args.year2}_{mass}.png")


        # Plot S, B vs BDT cut
        fig, ax1 = plt.subplots(figsize = (8,8))
        ax2 = ax1.twinx()
        title = f"S and B in SR vs BDT Cut for m_{{a}} = {mass.replace('_',' ')}"
        ln1 = ax1.scatter(cut1+invalid_cut1, signal_sr1, s=5, c="b", label=f"{args.year1} S")
        ln2 = ax2.scatter(cut1+invalid_cut1, bkg_sr1, s=5, c="r", label=f"{args.year1} B")
        vln1 = ax1.axvline(max_cut1, color="g", ls=":", label=f"{args.year1} Optimal Cut")
        ln3 = ax1.scatter(cut2+invalid_cut2, signal_sr2, s=5, c="tab:cyan", marker="s", label=f"{args.year2} S")
        ln4 = ax2.scatter(cut2+invalid_cut2, bkg_sr2, s=5, c="tab:orange", marker="s", label=f"{args.year2} B")
        vln2 = ax1.axvline(max_cut2, color="k", ls="--", label=f"{args.year2} Optimal Cut")
        ax1.set_xlabel("Cut on BDT Score", fontsize=14)
        ax1.set_ylabel("Signal in Signal Region (SR)", fontsize=14)
        ax2.set_ylabel("Background in Signal Region (SR)", fontsize=14)
        ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
        ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax2.xaxis.set_minor_locator(AutoMinorLocator(10))
        all_lns = [ln1, ln2, vln1, ln3, ln4, vln2]
        ax1.legend(all_lns, [l.get_label() for l in all_lns], fontsize=14, loc=9)
        ax1.set_ylim(top=198.0, bottom=28.0)
        
        # Add text to plot
        ax1.text(0.465, 140.0, "BDT Cut", weight="bold", fontsize=13)
        ax1.text(0.435, 135.0, f"{args.year1}: {max_cut1:.4f}", fontsize=13)
        ax1.text(0.435, 130.0, f"{args.year2}: {max_cut2:.4f}", fontsize=13)
        ax1.text(0.43, 120.0, "Optimal AMS", weight="bold", fontsize=13)
        ax1.text(0.45, 115.0, f"{args.year1}: {max_ams1:.2f}", fontsize=13)
        ax1.text(0.45, 110.0, f"{args.year2}: {max_ams2:.2f}", fontsize=13)
        ax1.text(0.35, 100.0, "Signal at Optimal AMS", weight="bold", fontsize=13)
        ax1.text(0.45, 95.0, f"{args.year1}: {signal_sr1[len(cut1)-1]:.2f}", fontsize=13)
        ax1.text(0.45, 90.0, f"{args.year2}: {signal_sr2[len(cut2)-1]:.2f}", fontsize=13)
        ax1.text(0.38, 80.0, "Bkg at Optimal AMS", weight="bold", fontsize=13)
        ax1.text(0.45, 75.0, f"{args.year1}: {bkg_sr1[len(cut1)-1]:.2f}", fontsize=13)
        ax1.text(0.45, 70.0, f"{args.year2}: {bkg_sr2[len(cut2)-1]:.2f}", fontsize=13)

        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.tight_layout()
        if not os.path.exists(f"{cwd}/plots/BDT/ams_vs_bdt_{args.year1}_{args.year2}/"):
            Path(f"{cwd}/plots/BDT/ams_vs_bdt_{args.year1}_{args.year2}/").mkdir(parents=True, exist_ok=True)
        plt.savefig(f"{cwd}/plots/BDT/ams_vs_bdt_{args.year1}_{args.year2}/sb_vs_bdt_{args.year1}_{args.year2}_{mass}.png")
        print(f"Saved plot sb_vs_bdt_{args.year1}_{args.year2}/ams_vs_bdt_{args.year1}_{args.year2}_{mass}.png")



if __name__ == "__main__":
    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser", conflict_handler="resolve")
    parser.add_argument("-i1", "--input1", required=True, help="Load first file for AMS vs BDT.", type=str)
    parser.add_argument("-i2", "--input2", required=True, help="Load second file for AMS vs BDT.", type=str)
    parser.add_argument("-y1", "--year1", required=True, help="Year of first sample.", type=str)
    parser.add_argument("-y2", "--year2", required=True, help="Year of second sample.", type=str)
    args = parser.parse_args()

    compareAMS(args)
