#!/usr/bin/env python3
 
import ROOT
import awkward as ak
import numpy as np
import json
 
ROOT.gROOT.SetBatch(True)
 
def getEffSigma(_h):
 
    nbins = _h.GetXaxis().GetNbins()
    binw = _h.GetXaxis().GetBinWidth(1)
    xmin = _h.GetXaxis().GetXmin()
 
    mu = _h.GetMean()
    rms = _h.GetRMS()
    total = _h.Integral()
 
    nWindow = (
        int(rms / binw)
        if (rms / binw) < 0.1 * nbins
        else int(0.1 * nbins)
    )
 
    rlim = 0.683 * total
    wmin = 9999999.
 
    for iscan in range(-nWindow, nWindow + 1):
 
        i_centre = int((mu - xmin) / binw + 1 + iscan)
 
        if i_centre < 1 or i_centre > nbins:
            continue
 
        x_centre = (i_centre - 0.5) * binw + xmin
 
        x_up = x_centre
        x_down = x_centre
 
        i_up = i_centre
        i_down = i_centre
 
        y = _h.GetBinContent(i_centre)
        r = y
 
        reachedLimit = False
 
        for j in range(1, nbins):
 
            if reachedLimit:
                continue
 
            if (i_up < nbins) and (not reachedLimit):
 
                i_up += 1
                x_up += binw
 
                y = _h.GetBinContent(i_up)
                r += y
 
                if r > rlim:
                    reachedLimit = True
 
            else:
                print(f"Reached upper edge for {_h.GetName()}")
                return 0
 
            if not reachedLimit:
 
                if i_down > 1:
 
                    i_down -= 1
                    x_down -= binw
 
                    y = _h.GetBinContent(i_down)
                    r += y
 
                    if r > rlim:
                        reachedLimit = True
 
                else:
                    print(f"Reached lower edge for {_h.GetName()}")
                    return 0
 
        if y == 0:
            dx = 0.
        else:
            dx = (r - rlim) * (binw / y)
 
        w = (x_up - x_down + binw - dx) * 0.5
 
        if w < wmin:
            wmin = w
 
    return wmin
 
def compute_sigma_from_pdf(pdf, mass):
 
    hist = pdf.createHistogram(
        "h_sig",
        mass,
        ROOT.RooFit.Binning(3200)
    )
 
    hist.Scale(80.0 / 160.0)
 
    return getEffSigma(hist)
 
def compute_peak_and_hwhm(pdf, mass):
 
    scan_min = 115.0
    scan_max = 135.0
    scan_step = 0.01
 
    masses = []
    values = []
 
    obs = ROOT.RooArgSet(mass)
 
    m = scan_min
 
    while m <= scan_max + 1e-9:
 
        mass.setVal(m)
 
        value = pdf.getVal(obs)
 
        masses.append(m)
        values.append(value)
 
        m += scan_step
 
    masses = np.asarray(masses)
    values = np.asarray(values)
 
    #finding the peak value of signal mdels
    peak_index = np.argmax(values)
 
    peak_mass = masses[peak_index]
    peak_value = values[peak_index]
 
    hist = pdf.createHistogram(
        "h_fwhm",
        mass,
        ROOT.RooFit.Binning(3200)
    )
 
    half_max = 0.5 * hist.GetMaximum()
 
    low_bin = hist.FindFirstBinAbove(half_max)
    high_bin = hist.FindLastBinAbove(half_max)
 
    fwhm_low = hist.GetBinCenter(low_bin)
    fwhm_high = hist.GetBinCenter(high_bin)
 
    fwhm = fwhm_high - fwhm_low
    hwhm = fwhm / 2.0
 
    return peak_mass, peak_value, fwhm, hwhm
 
 
def process_workspace(file_path, pdf_name):
 
    f = ROOT.TFile.Open(file_path)
 
    ws = f.Get("wsig_13TeV")
 
    pdf = ws.pdf(pdf_name)
 
    mass = ws.var("CMS_hgg_mass")
 
    mh = ws.var("MH")
    mh.setVal(125.0)
 
    sigma_eff = compute_sigma_from_pdf(pdf, mass)
 
    peak_mass, peak_value, fwhm, hwhm = compute_peak_and_hwhm(
        pdf,
        mass
    )
 
    f.Close()
 
    return {
        "sigma_eff": round(sigma_eff, 2),
        "peak_mass": round(peak_mass, 2),
        "peak_value": round(peak_value, 2),
        "fwhm": round(fwhm, 2) if fwhm is not None else None,
        "hwhm": round(hwhm, 2) if hwhm is not None else None
    }
 
 
results = {
    "Run3_2024": {},
    "Run2": {
        "2016": {},
        "2017": {},
        "2018": {}
    }
}
 
 
for m in range(15, 62, 5):
 
    file_path = f"/eos/home-t/takumar/gitcode/higgsdna_finalfits_tutorial_24/07_FinalFits/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/outdir_bdt-{m}-2024_lcg3/signalFit/output/CMS-HGG_sigfit_bdt-{m}-2024_lcg3_ma{m}_2024_CAT1.root"
 
    pdf_name = "hggpdfsmrel_GG2H_2024_CAT1_13TeV"
 
    print(f"Processing Run 3, mass = {m} GeV")
 
    results["Run3_2024"][str(m)] = process_workspace(
        file_path,
        pdf_name
    )
 
 
for year in [2016, 2017, 2018]:
 
    for m in range(15, 62, 5):
 
        file_path = (
            f"/eos/home-t/takumar/gitcode/"
            f"higgsdna_finalfits_tutorial_24/07_FinalFits/"
            f"CMSSW_14_1_0_pre4/src/flashggFinalFit/"
            f"Combine/badder/datacards/Datacards/Models/"
            f"{year}/"
            f"CMS-HGG_mva_13TeV_M{m}_sigfit.root"
        )
 
        pdf_name = f"hggpdfsmrel_{year}_13TeV_H4G_Cat0"
 
        print(f"Processing Run 2, year = {year}, mass = {m} GeV")
 
        results["Run2"][str(year)][str(m)] = process_workspace(
            file_path,
            pdf_name
        )
 
 
output_file = "sigma_results.json" #results that can be used later
 
with open(output_file, "w") as f:
    json.dump(results, f, indent=4)
 
print()
print("=" * 60)
print(f"Results saved to: {output_file}")
print("=" * 60)
