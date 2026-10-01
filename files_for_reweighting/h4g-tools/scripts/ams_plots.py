import os
import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

if __name__ == "__main__":
    
    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser", conflict_handler="resolve")
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    parser.add_argument("-xgb", "--xgb_name", type=str, required=True, help="XGBoost model name.")
    parser.add_argument("-zero", "--includeZero", action="store_true", help="Include zeroeth category on plot.")
    args = parser.parse_args()

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    cats = {}
    # Read in category information from all files
    for nCats in range(len([f for f in os.listdir(os.path.join(cwd, "cats", args.year, args.xgb_name)) if ("cats" in f and ".txt" in f and "bdt" not in f)])):
        with open(os.path.join(cwd, "cats", args.year, args.xgb_name, f"cats{nCats+1}.txt"), "r") as f:
            cats.update({nCats+1: f.readlines()})

    ams = {}
    # Store AMS for each mass point for each # of categories
    for key in cats.keys():
        for idx, mass in enumerate([f"{m} GeV" for m in range(15,65,5)]):
            if mass not in ams.keys():
                ams.update({mass: {}})

            # Get line with AMS in it without wasting compute on a loop
            ams_line = 1 + (key + 3) * idx
            line  = cats[key][ams_line]

            assert mass in line, (mass, line)
            assert f"# of cats: {key}" in line, (mass, line)

            ams[mass].update({key: float(line[line.find("AMS")+4: line.find("Signal Total")])})

    for idx, mass in enumerate(ams.keys()):
        x = [0] if args.includeZero else []
        x.extend(ams[mass].keys())

        y = [0.0] if args.includeZero else []
        y.extend(ams[mass].values())

        fig = plt.figure(idx)
        plt.plot(x,y,"-o", color="b")

        ax = fig.get_axes()[0]
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f' if "2024" in args.year else '%.5f'))
        plt.yticks(rotation=45)

        plt.xlabel("Number of categories", fontsize=12)
        plt.ylabel("AMS", fontsize=12)
        plt.grid()

        plt.locator_params(axis="x", integer=True)
        plt.locator_params(axis="y", integer=False)

        # Save figure
        savePath = f"{cwd}/plots/BDT/{args.year}/cats/{args.xgb_name}/{mass.replace(' ','_')}_categories.png"
        plt.savefig(savePath, bbox_inches='tight')
        print(f"Saved AMS plot for {mass}.")
        assert os.path.isfile(savePath), "Saved plot does not exist."
