import os
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

if __name__ == "__main__":
    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    mean = {
        "gauss1": [124.4219, 123.9407, 124.2483, 124.6320, 124.8286, 124.2900, 124.2565, 125.2223, 124.5921, 124.3193],
        "gauss2": [124.5735, 125.0381, 124.8130, 123.9007, 123.5775, 124.7303, 124.6999, 123.6266, 124.2935, 124.5953],
        "gauss3": [123.7177, 124.1710, 124.1832, 125.2875, 124.4438, 124.1844, 124.1751, 123.9859, 123.9579, 123.5108]
    }

    sigma = {
        "gauss1": [2.3281, 1.9039, 2.3228, 1.6772, 1.5493, 2.2326, 2.3481, 1.7163, 1.6208, 2.4613],
        "gauss2": [1.3921, 1.6327, 1.2968, 3.3012, 2.2838, 1.4990, 1.3911, 2.0464, 2.6610, 1.5152],
        "gauss3": [4.3208, 4.1559, 4.2226, 5.4147, 4.0397, 4.1993, 4.0842, 4.5768, 3.9635, 4.2174]
    }

    frac = {
        "gauss1": [0.4771, 0.3649, 0.4619, 0.5395, 0.4240, 0.3908, 0.4287, 0.4032, 0.4360, 0.4271],
        "gauss2": [0.2592, 0.2937, 0.2614, 0.3842, 0.2413, 0.3015, 0.2696, 0.3564, 0.2724, 0.3092],
        "gauss3": [0.2637, 0.3414, 0.2767, 0.0763, 0.3346, 0.3077, 0.3017, 0.2404, 0.2916, 0.2638]
    }

    # MEAN
    fig1, ax1 = plt.subplots(figsize = (6,6))
    x = list(range(15,65,5))
    plt.plot(x, mean["gauss1"], "-o", color="r")
    plt.plot(x, mean["gauss2"], "-o", color="b")
    plt.plot(x, mean["gauss3"], "-o", color="g")
    plt.xlabel("m(a) [GeV]", fontsize=14)
    plt.ylabel("Gaussian Mean ($\mu$) [GeV]", fontsize=14)

    ax = plt.gca()
    ymin = 123.0
    ymax = 126.5
    ax.set_ylim(ymin, ymax)
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.legend(labels=["Gaussian 1", "Gaussian 2", "Gaussian 3"], frameon=True, loc="upper right", fontsize=14)
    plt.savefig(f"{cwd}/plots/mean_summary-sameN.pdf")
    print(f"Saved plot to {cwd}/plots/mean_summary-sameN.pdf")

    # WIDTH
    fig2, ax2 = plt.subplots(figsize = (6,6))
    x = list(range(15,65,5))
    plt.plot(x, sigma["gauss1"], "-o", color="r")
    plt.plot(x, sigma["gauss2"], "-o", color="b")
    plt.plot(x, sigma["gauss3"], "-o", color="g")
    plt.xlabel("m(a) [GeV]", fontsize=14)
    plt.ylabel("Gaussian Width ($\sigma$) [GeV]", fontsize=14)

    ax = plt.gca()
    ymin = 1.0
    ymax = 6.5
    ax.set_ylim(ymin, ymax)
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.legend(labels=["Gaussian 1", "Gaussian 2", "Gaussian 3"], frameon=True, loc="upper right", fontsize=14)
    plt.savefig(f"{cwd}/plots/sigma_summary-sameN.pdf")
    print(f"Saved plot to {cwd}/plots/sigma_summary-sameN.pdf")

    fig3, ax3 = plt.subplots(figsize = (6,6))
    x = list(range(15,65,5))
    plt.plot(x, frac["gauss1"], "-o", color="r")
    plt.plot(x, frac["gauss2"], "-o", color="b")
    plt.plot(x, frac["gauss3"], "-o", color="g")
    plt.xlabel("m(a) [GeV]", fontsize=14)
    plt.ylabel("Fractional Contribution", fontsize=14)

    ax = plt.gca()
    ymin = 0.0
    ymax = 0.75
    ax.set_ylim(ymin, ymax)
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.legend(labels=["Gaussian 1", "Gaussian 2", "Gaussian 3"], frameon=True, loc="upper right", fontsize=14)
    plt.savefig(f"{cwd}/plots/fraction_summary-sameN.pdf")
    print(f"Saved plot to {cwd}/plots/fraction_summary-sameN.pdf")

