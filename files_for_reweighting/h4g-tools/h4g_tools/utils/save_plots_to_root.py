import os
from ROOT import TFile, TH1D

def savePlots(
    directory: str,
    hists: dict,
    names: list,
    save_path: str
) -> None:
    """
    Wrapper to save returned histograms to a ROOT file for use in macros.
    """

    copy_hists = []

    # Update name only - let title be the name for the legends (when needed)
    double_loop = False
    if "GeV" in list(hists.keys())[0]:
        if type(hists[list(hists.keys())[0]]) is not TH1D:
            double_loop = True
    if not double_loop:
        for h, name in zip(hists.values(), names):
            h.SetName(name)
            copy_hists.append(h)
    elif double_loop:
        for m, d in hists.items():
            for h, name in zip(d.values(), names):
                new_name = f"{name}_{m.replace('Signal_', '')}"
                h.SetName(new_name)
                copy_hists.append(h)

    f = TFile.Open(save_path, "UPDATE")
    d = f.GetDirectory(directory)
    if not d:
        d = f.mkdir(directory)
    d.cd()
    for h in copy_hists:
        h.Write()

    f.Close()

