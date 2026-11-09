import xgboost as xgb
import matplotlib.pyplot as plt

def inspect_xgb_model(model_file_path, output_prefix="tree_plot"):
    print(f"=== Inspecting: {model_file_path} ===")
    bst = xgb.Booster()
    bst.load_model(model_file_path)

    # Total number of trees
    trees_dump = bst.get_dump()
    num_trees = len(trees_dump)
    print(f"n_estimators (Total Trees): {num_trees}")

    # Model parameters
    config = bst.save_config()
    print("\nModel Configuration / Parameters:")
    print(config)

    # Tree structure preview
    print("\nTree Structure Preview (First 2 Trees):")
    for i in range(min(2, num_trees)):
        print(f"--- Tree {i} ---")
        print(trees_dump[i])

    # Save tree diagrams to PNG (e.g., saving Tree 0)
    if num_trees > 0:
        fig, ax = plt.subplots(figsize=(20, 10))
        # num_trees=0 plots the first tree (index 0)
        xgb.plot_tree(bst, num_trees=0, ax=ax)
        
        output_png = f"{output_prefix}_tree0.png"
        plt.savefig(output_png, bbox_inches='tight', dpi=300)
        plt.close(fig)
        print(f"\nSaved tree 0 diagram to: {output_png}\n")

# Inspect and save diagrams for both models
inspect_xgb_model(
    '/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/files_for_reweighting/h4g-tools/models/ArcTrialBDT_12Oct2024.xgb',
    output_prefix="ArcTrialBDT_12Oct2024"
)

inspect_xgb_model(
    '/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/files_for_reweighting/h4g-tools/models/ArcTrialBDTwInputParams_13Oct2024.xgb',
    output_prefix="ArcTrialBDTwInputParams_13Oct2024"
)
