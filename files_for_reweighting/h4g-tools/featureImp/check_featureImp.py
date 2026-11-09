import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import xgboost as xgb

model = xgb.XGBClassifier()
model.load_model("/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/files_for_reweighting/h4g-tools/models/ArcTrialBDT_12Oct2024.xgb")

fig, ax = plt.subplots(figsize=(10, 8))
xgb.plot_importance(
    model, importance_type="weight", max_num_features=10, ax=ax
)

plt.tight_layout()
plt.savefig("feature_importance.png", dpi=300, bbox_inches="tight")
plt.close()
