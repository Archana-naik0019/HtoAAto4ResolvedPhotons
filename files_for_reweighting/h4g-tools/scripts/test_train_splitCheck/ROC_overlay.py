import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

models = {
    "50/50": "/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/files_for_reweighting/h4g-tools/scripts/plots/BDT/2024/train_checks/ROC_curves_2024.npz",
    "80/20": "/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/files_for_reweighting/h4g-tools/scripts/plots/BDT/2024/train80_20_checks/ROC_curves_2024.npz",
}
outdir = os.path.join(os.path.dirname(models["80/20"]), "split_comparison")
os.makedirs(outdir, exist_ok=True)

# Subdirectories for unzoomed plots
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

# one colour per model, solid = test, dashed = train
col = {"50/50": "tab:blue", "80/20": "tab:red"}
ls = {"test": "-", "train": "--"}

def overlay(m, fname, title, log_scale=True):
    fig, ax = plt.subplots(figsize=(7, 6.5))
    for k in models:
        for split in ["train", "test"]:
            if not has(k, split, m):
                continue
            f, t = curve(k, split, m)
            ax.plot(f, t, color=col[k], ls=ls[split], lw=1.5,
                    label=f"{k} {split} (AUC={auc(k, split, m):.5f})")
    
    if log_scale:
        style_log(ax)
        save_path = f"{outdir}/{fname}"
    else:
        style_linear(ax)
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

# --- Generate Linear Scale (Unzoomed) Plots ---
overlay("net", "ROC_net_split_comparison_unzoomed", "Net ROC (Unzoomed): 50/50 vs 80/20 split", log_scale=False)

pm_unzoomed = f"{unzoomed_outdir}/per_mass"; os.makedirs(pm_unzoomed, exist_ok=True)
for m in masses:
    overlay(m, f"per_mass/ROC_{m}GeV_split_comparison_unzoomed", f"{m} GeV (Unzoomed): 50/50 vs 80/20", log_scale=False)

# AUC summary table: train-test gap per model
print(f"{'mass':>6} | " + " | ".join(f"{k} train  test   gap" for k in models))
for m in masses + ["net"]:
    cells = []
    for k in models:
        if has(k, "train", m) and has(k, "test", m):
            a, b = auc(k, "train", m), auc(k, "test", m)
            cells.append(f"{a:.5f} {b:.5f} {a-b:+.5f}")
        else:
            cells.append("n/a")
    print(f"{str(m):>6} | " + " | ".join(cells))
