import matplotlib.pyplot as plt
import os

if __name__ == "__main__":

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Plot HLT bit efficiency from signal MC
    fig, ax1 = plt.subplots(figsize = (6,6))
    title = f"Ratio of Selection Efficiencies"
    x = list(range(15,65,5))
    y = [1.043,
         1.050,
         1.054,
         1.056,
         1.060,
         1.065,
         1.069,
         1.069,
         1.066,
         1.066,
    ]
    plt.plot(x, y, "-o")
    plt.xlabel("m(a) [GeV]", fontsize=14)
    plt.ylabel("$Eff_{no\ HLT\ bit}$ / $Eff_{with\ HLT\ bit}$", fontsize=14)
    plt.grid()

    ax = plt.gca()
    ymin = 1.0
    ymax = 1.1
    ax.set_ylim(ymin, ymax)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig(f"{cwd}/plots/standard/2024/selection_eff_with_HLT.pdf")
    print(f"Saved plot to {cwd}/plots/standard/2024/selection_eff_with_HLT.pdf")

