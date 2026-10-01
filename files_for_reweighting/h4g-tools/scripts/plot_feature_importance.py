import os
import numpy as np
import matplotlib.pyplot as plt
import xgboost as xgb
from matplotlib.ticker import AutoMinorLocator

if __name__ == "__main__":
    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    #####################################################
    ##### !!! Edit these paths for your setup !!! #######
    model_path = f"{cwd}/../models/ES_BDT_2024.xgb"
    plot_path = f"{cwd}/plots/BDT/2024/"
    #####################################################

    features = (
        "$a_1$ $p_T$",
        "$a_2$ $p_T$",
        "$\\frac{\Delta R_{aa}}{m_{4\gamma}}$",
        "$\\frac{m_{a1} - m_{Hyp}}{m_{4\gamma}}$",
        "$\\frac{m_{a2} - m_{Hyp}}{m_{4\gamma}}$",
        "$m_{a1} - m_{a2}$",
        "$cos(\\theta_{a\gamma})$",
        "$\gamma_1$  mvaID",
        "$\gamma_2$  mvaID",
        "$\gamma_3$  mvaID",
        "$\gamma_4$  mvaID",
    )

    fimp_keys = [f"f{i}" for i in range(len(features))]
    fimp_map = {features[i]: fimp_keys[j] for i,j in zip(range(len(features)), [4,5,6,7,8,9,10,0,1,2,3])}

    bdt = xgb.XGBClassifier()
    bdt.load_model(model_path)
    booster = bdt.get_booster()
    y_pos = np.arange(len(features))

    for imp in ["gain", "cover", "weight"]:
        fig, ax = plt.subplots()
        fimp = booster.get_score(importance_type=imp)
        values = [fimp[fimp_map[k]] for k in features]

        ax.barh(y_pos, values, align='center')
        ax.set_yticks(y_pos, labels=features)
        ax.set_xlabel(f"Feature Importance ({imp})")
        ax.invert_yaxis()
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        plt.tight_layout()

        fig.savefig(f"{plot_path}feature_importance_{imp}.pdf")
        print(f"Saved plot to {plot_path}/feature_importance_{imp}.pdf")


