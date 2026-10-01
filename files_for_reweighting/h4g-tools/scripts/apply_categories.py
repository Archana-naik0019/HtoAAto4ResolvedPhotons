import os
import json
from h4g_tools.utils.loading import pathSetup, loadPerMass
from h4g_tools.utils.runner_utils import get_parser_cats


def split_twos(
    lst: list
) -> list:
    
    new_list = []
    for idx in range(0,len(lst), 2):
        new_list.append([lst[idx], lst[idx+1]])

    return new_list


if __name__ == "__main__":
    
    # Run pathSetup and get parser
    parser = get_parser_cats()
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    parser.add_argument("-xgb", "--xgb_name", type=str, required=True, help="XGBoost model name.")
    parser.add_argument("-bkg", "--bkg", action="store_true", required=False, help="Create BDT efficiency JSON for background samples.")
    args = parser.parse_args()

    pathSetup(year=args.year, model=args.xgb_name)

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    if not args.bkg:
        print("\nLoading signal samples:")
        samples = loadPerMass(f"{cwd}/../../scripts/outputs/BDT/{args.year}/{args.xgb_name}/signal/", ["BDT_score"], bdt=True)
        print()
    elif args.bkg:
        print("\nLoading background samples:")
        samples = loadPerMass(f"{cwd}/../../scripts/outputs/BDT/{args.year}/{args.xgb_name}/background/", ["BDT_score"], bdt=True)
        print()

    nCats = args.nCats
    cats = dict.fromkeys(samples.keys(), list().copy())

    with open(f"{cwd}/cats/{args.year}/{args.xgb_name}/cats{nCats}.txt", "r") as f:
        mass_pt = None
        nCats_file = -1
        for line in f:
            if "Mass Point: " in line:
                mass_point = line[line.find("GeV")-3 : line.find("GeV")+3].strip().replace(" ","_")
                nCats_file = int(line[line.find("# of cats:")+10: line.find("# of cats:")+12])

            if "Cuts: " in line and nCats_file == nCats:
                cuts_str = line[len("Cuts: "):]
                cats[mass_point] = sorted([float(c) for c in cuts_str.split(",")], reverse=True)
                assert len(cats[mass_point]) % 2 == 0, len(cats[mass_point])

                cats[mass_point] = split_twos(cats[mass_point])

    assert all([len(l) == nCats for l in cats.values()]), [len(l) for l in cats.values()]

    new_samples = {}
    events = {}
    for mass in [f"{m}_GeV" for m in range(15,65,5)]:
        print(f"~~~  Mass Point: {mass.replace('_',' ')}  ~~~")
        new_samples.update({mass: {}})
        events.update({mass: {}})
        totals = len(samples[mass])
        for i in range(len(cats[mass])):
            if i not in new_samples[mass].keys():
                new_samples[mass].update({i: None})
            new_samples[mass][i] = samples[mass][(samples[mass]["BDT_score"] <= cats[mass][i][0]) & (samples[mass]["BDT_score"] > cats[mass][i][1])]
            events[mass].update({i: len(new_samples[mass][i])})
            events[mass].update({f"cut{i}": (cats[mass][i][1], cats[mass][i][0])})

            if i == nCats:
                break
            
            eff = len(new_samples[mass][i])/totals * 100
            eff = f"{eff: .2f}"
            print(f"Events after cut for Cat {i+1} {tuple(cats[mass][i])}: {len(new_samples[mass][i])} of {totals} ({eff.strip()}%)")
        print()

    path_out = f"{cwd}/cats/{args.year}/{args.xgb_name}/bdt_{args.year}_cats{args.nCats}.json" 
    path_out = path_out if not args.bkg else path_out.replace("_cats", "bkg_cats")
    with open(path_out, "w") as f:
        json.dump(events, f, indent=4)

    print(f"\nSaved dictionary to {path_out}")

