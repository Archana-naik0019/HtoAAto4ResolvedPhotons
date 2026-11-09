import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Current working directory as output path
outdir = os.getcwd()

models = {
    "50/50": "/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/files_for_reweighting/h4g-tools/scripts/plots/BDT/2024/train_checks/ROC_curves_2024.npz",
    "80/20": "/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/files_for_reweighting/h4g-tools/scripts/plots/BDT/2024/train_newPresel8020/ROC_curves_2024.npz",
}

# Subdirectories created inside the current directory
unzoomed_outdir = os.path.join(outdir, "unzoomed")
os.makedirs(unzoomed_outdir, exist_ok=True)

data = {k: np.load(v) for k, v in models.items()}
masses = list(range(15, 65, 5))
trapz = getattr(np, "trapezoid", None) or np.trapz

def has(k, split, m):
    return f"{split}_{m}_fpr" in data[k].files

def curve(k, split, m):
    return data[k][f"{split}_{m}_fpr"], data[k][f"{split}_{m}_tpr"]

def auc(k, split, m):
    f, t = curve(k, split, m)
    return trapz(t, f)

def style_log(ax, ymin=0.65):
    x = np.logspace(-4, 0, 200)
    ax.plot(x, x, "k:", lw=1, label="Random guess")
    ax.set_xscale("log"); ax.set_xlim(1e-4, 1); ax.set_ylim(ymin, 1.005)
    ax.set_xlabel("Background efficiency (log scale)")
    ax.set_ylabel("Signal efficiency")
    ax.grid(True, which="both", alpha=0.3)

def style_linear(ax, pad=0.02):
    x = np.linspace(0, 1, 200)
    ax.plot(x, x, "k:", lw=1, label="Random guess")
    ax.set_xlim(0 - pad, 1 + pad)
    ax.set_ylim(0 - pad, 1 + pad)
    ax.set_xlabel("Background efficiency (FPR)")
    ax.set_ylabel("Signal efficiency (TPR)")
    ax.grid(True, which="both", alpha=0.3)

# Style for 4-curve overlaid plots: 80/20 = solid line (-), 50/50 = dashed line (--)
col = {"50/50": "tab:blue", "80/20": "tab:red"}
ls = {"50/50": "--", "80/20": "-"}

# Custom style for single-model individual plots: Train=Blue (solid), Test=Red (solid)
col_single = {"train": "tab:blue", "test": "tab:red"}
ls_single = {"train": "-", "test": "-"}

def overlay(m, fname, title, log_scale=True, selected_models=None, single_model_style=False):
    if selected_models is None:
        selected_models = list(models.keys())

    fig, ax = plt.subplots(figsize=(7, 6.5))
    for k in selected_models:
        for split in ["train", "test"]:
            if not has(k, split, m):
                continue
            f, t = curve(k, split, m)
            
            # Apply formatting rules
            if single_model_style:
                c = col_single[split]
                l_style = ls_single[split]
            else:
                c = col[k]
                l_style = ls[k]  # Solid for 80/20, Dashed for 50/50
            
            ax.plot(f, t, color=c, ls=l_style, lw=1.5,
                    label=f"{k} {split} (AUC={auc(k, split, m):.5f})")

    if log_scale:
        style_log(ax)
        save_path = f"{outdir}/{fname}"
    else:
        style_linear(ax, pad=0.02)
        save_path = f"{unzoomed_outdir}/{fname}"

    ax.set_title(title)
    ax.legend(fontsize=8, loc="lower right")
    fig.savefig(f"{save_path}.pdf"); fig.savefig(f"{save_path}.png", dpi=150)
    plt.close(fig)

# --- Generate Log Scale (Zoomed) Plots ---
overlay("net", "ROC_net_split_comparison", "Net ROC: 50/50 vs 80/20 split (mixed bkg)", log_scale=True)

pm = f"{outdir}/per_mass"; os.makedirs(pm, exist_ok=True)
for m in masses:
    overlay(m, f"per_mass/ROC_{m}GeV_split_comparison", f"{m} GeV: 50/50 vs 80/20 (mixed bkg)", log_scale=True)

# --- Generate Linear Scale Plots ---
overlay("net", "ROC_net_split_comparison_unzoomed", "Net ROC: 50/50 vs 80/20 split", log_scale=False)

pm_unzoomed = f"{unzoomed_outdir}/per_mass"; os.makedirs(pm_unzoomed, exist_ok=True)
for m in masses:
    overlay(m, f"per_mass/ROC_{m}GeV_split_comparison_unzoomed", f"{m} GeV: 50/50 vs 80/20", log_scale=False)

# --- Individual Single-Model Linear Net Plots ---
overlay("net", "ROC_net_50_50_unzoomed", "Net ROC: 50/50 Train vs Test", 
        log_scale=False, selected_models=["50/50"], single_model_style=True)

overlay("net", "ROC_net_80_20_unzoomed", "Net ROC: 80/20 Train vs Test", 
        log_scale=False, selected_models=["80/20"], single_model_style=True)

# --- Individual Per-Mass Single-Model Linear Plots ---
for m in masses:
    overlay(m, f"per_mass/ROC_{m}GeV_50_50_unzoomed", f"{m} GeV: 50/50 Train vs Test", 
            log_scale=False, selected_models=["50/50"], single_model_style=True)
    overlay(m, f"per_mass/ROC_{m}GeV_80_20_unzoomed", f"{m} GeV: 80/20 Train vs Test", 
            log_scale=False, selected_models=["80/20"], single_model_style=True)

# AUC summary table: train-test gap per model
print(f"{'mass':>6} | " + " | ".join(f"{k} train   test    gap" for k in models))
for m in masses + ["net"]:
    cells = []
    for k in models:
        if has(k, "train", m) and has(k, "test", m):
            a, b = auc(k, "train", m), auc(k, "test", m)
            cells.append(f"{a:.5f} {b:.5f} {a-b:+.5f}")
        else:
            cells.append("n/a")
    print(f"{str(m):>6} | " + " | ".join(cells))
