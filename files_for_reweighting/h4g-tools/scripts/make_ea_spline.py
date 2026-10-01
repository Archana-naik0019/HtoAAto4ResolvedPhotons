import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

if __name__ == "__main__":
    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    numer = np.array([35133, 34716, 34252, 34164, 38151, 38043, 43989, 59189, 79816, 98538])
    denom = np.array([1025691, 971967, 969906, 981050, 1047791, 981828, 974661, 1022301, 994105, 955078])
    y = numer / denom
    assert len(numer) == len(denom), (len(numer), len(denom))
    
    spline = np.zeros(63-15)  # Need 62+1 becuase of numpy being stupid with indices
    for idx in range(len(spline)):
        if (15 + idx) % 5 == 0:
            if 15 + idx < 60:
                m = (y[int((idx/5)+1)] - y[int(idx/5)]) / 5
                b = y[int(idx/5)]
                spline[idx] = y[int(idx/5)]
            else:
                b = y[int(idx/5)]
                spline[idx] = y[int(idx/5)]
            #print(15 + idx, y[int(idx/5)], (m, idx, b), end="\t")
            print(15 + idx, end="\t")
        else:
            spline[idx] = m*(idx%5) + b
        print(spline[idx])

    fig1, ax1 = plt.subplots(figsize = (6,6))
    masses = list(range(15,63,1))
    plt.plot(masses, spline, "-o", color="b", ms=3.6)
    plt.plot(list(range(15,65,5)), y, "o", color="r", ms=3.6)
    plt.xlabel("m(a) [GeV]", fontsize=14)
    plt.ylabel("Event Selection Efficiency", fontsize=14)

    ax = plt.gca()
    ymin = 0.025
    ymax = 0.135
    ax.set_ylim(ymin, ymax)
    ax.yaxis.set_major_locator(MultipleLocator(0.01))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.legend(labels=["Interpolated", "Nominal"], frameon=True, loc="upper right", fontsize=14)
    plt.savefig(f"{cwd}/plots/ea_spline_overlay.pdf")
    print(f"Saved plot to {cwd}/plots/ea_spline_overlay.pdf")

