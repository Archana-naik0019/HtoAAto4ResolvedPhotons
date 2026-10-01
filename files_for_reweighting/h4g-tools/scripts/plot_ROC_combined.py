import json
import matplotlib.pyplot as plt
import os
import numpy as np
from typing import Optional


if __name__ == "__main__":

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    path_2022merged = f"{cwd}/plots/BDT/2022/ROC_2022.json"
    path_2022preEE = f"{cwd}/plots/BDT/2022preEE/ROC_2022preEE.json"
    path_2022postEE = f"{cwd}/plots/BDT/2022postEE/ROC_2022postEE.json"

    with open(path_2022merged, "r") as f:
        j2022merged = json.load(f)
    with open(path_2022preEE, "r") as f:
        j2022preEE = json.load(f)
    with open(path_2022postEE, "r") as f:
        j2022postEE = json.load(f)

    def make_plot(path: str, zoomed: Optional[bool] = False):
        plt.close()
        auc = []
        colors = ["red", "blue", "green"]
        for c, model, era in zip(colors, [j2022merged, j2022preEE, j2022postEE], ["2022 Merged", "2022 preEE", "2022 postEE"]):
            # FPR = false positive rate
            # TPR = true positive rate
            # Thresholds = cuts on BDT score
            # AUC = area under ROC curve
            # ROC = receiver operating characteristic curve
            fpr = np.array(model["train"]["fpr"])
            tpr = np.array(model["train"]["tpr"])
            auc.append(model["train"]["auc"])

            plt.plot(
                1-fpr,
                tpr,
                color = c,
                label = era,
                markersize = 1
            )

            print(era, model["train"]["events"]["signal"], model["train"]["events"]["bkg"])

        low = 0.0
        high = 0.2
        if not zoomed:
            plt.text(low + 0.71, high + 0.575 + 0.04, f"2022 Merged AUC: {auc[0]:.4}")
            plt.text(low + 0.71, high + 0.575 + 0.02, f"2022 preEE AUC: {auc[1]:.4}")
            plt.text(low + 0.71, high + 0.575, f"2022 postEE AUC: {auc[1]:.4}")

            plt.xlim(0.7, 1.05)
            plt.ylim(0.7, 1.05)
        elif zoomed:
            plt.xlim(0.9, 1.02)
            plt.ylim(0.9, 1.02)

        plt.xlabel("Signal Efficiency", fontsize=12)
        plt.ylabel("Background Rejection", fontsize=12)
        #plt.xlabel("Signal Efficiency (1-fpr)", fontsize=12)
        #plt.ylabel("Background Rejection (tpr)", fontsize=12)
        plt.legend(frameon=False)
        plt.title("Training ROC")

        # Save figure
        plt.savefig(path)
        assert os.path.isfile(path), "Saved plot does not exist."
        print(f"Saved {path.split('/')[-1]} to file.")


    path = f"{cwd}/plots/BDT/ROC_train_compare_eras.png"
    path_zoomed = f"{cwd}/plots/BDT/ROC_train_compare_eras_zoomed.png"
    make_plot(path)
    make_plot(path_zoomed, zoomed=True)
