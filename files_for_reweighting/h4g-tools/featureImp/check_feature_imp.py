import matplotlib

# Set non-interactive backend for headless environments like LXPlus
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import xgboost as xgb

# 1. Initialize and load model
model = xgb.XGBClassifier()
model.load_model("/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/files_for_reweighting/h4g-tools/models/train80_20_checks.xgb")

# 2. Plot feature importance
fig, ax = plt.subplots(figsize=(10, 8))
xgb.plot_importance(
    model, importance_type="weight", max_num_features=11, ax=ax
)

# 3. Adjust layout and save to PNG
plt.tight_layout()
plt.savefig("OldPresel_8020.png", dpi=300, bbox_inches="tight")
plt.close()

print("Feature importance plot successfully saved as feature_importance.png")
