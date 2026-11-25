"""
effSigma (sliding window, 68.3%) + peak from a RooFit sum-of-3-Gaussians, per signal mass point.

Peak procedure:
  1. TH1 of mass_gggg in [115, 135], 160 bins (weighted, Sumw2 on)
  2. Highest bin = first guess of the peak.
  3. RooFit fit of RooAddPdf(3 x RooGaussian) in [peak-4, peak+4] GeV.
       * fractions are recursive -> f1+f2+f3 = 1 by construction
       * sigma1 < sigma2 < sigma3 enforced via sigma2 = sigma1+d2, sigma3 = sigma2+d3 (d >= 0)
  4. Window is re-centred on the fitted peak. (*fitting is done based on NLL minimization using minuit)
  5. Reported peak = x where the fitted pdf is maximal (fine scan of the pdf, just additional).

Usage: python3 2sig_win_3gauss.py --apply-scaling --filter-mass --save-plots --no-info-box
"""

import argparse
import glob
import os

import awkward as ak
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

XLO, XHI = 115.0, 135.0
RF = ROOT.RooFit


# effective sigma: sliding window containing 68.3% of the (weighted) events
# ----------------------------------------------------------------------
def getEffSigma_unbinned(masses, weights=None, fraction=0.683):
    """Narrowest window holding `fraction` of the events. effSigma = (high - low)/2."""
    if weights is None:
        s = np.sort(masses)
        N = len(s)
        k = int(np.ceil(fraction * N))
        low_edges = s[: N - k + 1]
        high_edges = s[k - 1:]
        i = np.argmin(high_edges - low_edges)
        lo, hi = low_edges[i], high_edges[i]
    else:
        idx = np.argsort(masses)
        s = masses[idx]
        cw = np.cumsum(weights[idx])
        target = fraction * cw[-1]
        N = len(s)
        min_width, lo, hi = np.inf, 0.0, 0.0
        j = 0
        for i in range(N):
            prev = cw[i - 1] if i > 0 else 0.0
            while j < N and (cw[j] - prev) < target:
                j += 1
            if j >= N:
                break
            if s[j] - s[i] < min_width:
                min_width = s[j] - s[i]
                lo, hi = s[i], s[j]
    return 0.5 * (hi - lo), (lo, hi)


def effsigma_from_curve(xs, vals, fraction=0.683):
    """effSigma of a pdf sampled on a fine grid (narrowest interval with `fraction` of the area)."""
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (vals[1:] + vals[:-1]) * np.diff(xs))])
    cdf /= cdf[-1]
    j = np.searchsorted(cdf, cdf + fraction)
    ok = j < len(xs)
    widths = xs[j[ok]] - xs[ok]
    return 0.5 * widths.min()


def weight_diagnostics(masses, weights):
    n = len(masses)
    if weights is None:
        return dict(n=n, sumw=float(n), wmin=1.0, wmax=1.0, n_neg=0, neff=float(n))
    sw, sw2 = weights.sum(), (weights ** 2).sum()
    return dict(n=n, sumw=float(sw), wmin=float(weights.min()), wmax=float(weights.max()),
                n_neg=int((weights < 0).sum()), neff=float(sw ** 2 / sw2))


def scan_pdf(pdf, x, xs):
    obs = ROOT.RooArgSet(x)
    out = np.empty_like(xs)
    for i, xi in enumerate(xs):
        x.setVal(float(xi))
        out[i] = pdf.getVal(obs)
    return out


# one RooFit : 3 Gaussians, fractions sum to 1, several attempts
# ----------------------------------------------------------------------
def do_fit(h, center, eff_sigma, tag, half_win, n_tries, common_mean, seed):
    nbins = h.GetNbinsX()
    lo_fit, hi_fit = center - half_win, center + half_win

    x = ROOT.RooRealVar(f"x_{tag}", "m_{#gamma#gamma#gamma#gamma} [GeV]", XLO, XHI)
    x.setBins(nbins)
    x.setRange("fitRange", lo_fit, hi_fit)
    dh = ROOT.RooDataHist(f"dh_{tag}", "dh", ROOT.RooArgList(x), h)

    sg = max(eff_sigma, 0.3)
    mu1 = ROOT.RooRealVar(f"mu1_{tag}", "#mu_{1}", center, center - 1.5, center + 1.5)
    if common_mean:
        mu2 = mu3 = mu1
    else:
        mu2 = ROOT.RooRealVar(f"mu2_{tag}", "#mu_{2}", center, center - 2.5, center + 2.5)
        mu3 = ROOT.RooRealVar(f"mu3_{tag}", "#mu_{3}", center, center - 2.5, center + 2.5)
    # widths: sigma1 < sigma2 < sigma3 via positive increments; each capped at half_win/2 so that
    # sigma3 <= 1.5*half_win (a much wider Gaussian is flat inside the window = degenerate with the others)
    cap = 0.5 * half_win
    s1 = ROOT.RooRealVar(f"s1_{tag}", "#sigma_{1}", min(0.5 * sg, 0.9 * cap), 0.05, cap)
    d2 = ROOT.RooRealVar(f"d2_{tag}", "d_{2}", min(0.6 * sg, 0.9 * cap), 0.0, cap)
    d3 = ROOT.RooRealVar(f"d3_{tag}", "d_{3}", min(1.2 * sg, 0.9 * cap), 0.0, cap)
    s2 = ROOT.RooFormulaVar(f"s2_{tag}", "#sigma_{2}", "@0+@1", ROOT.RooArgList(s1, d2))
    s3 = ROOT.RooFormulaVar(f"s3_{tag}", "#sigma_{3}", "@0+@1+@2", ROOT.RooArgList(s1, d2, d3))

    g1 = ROOT.RooGaussian(f"g1_{tag}", "g1", x, mu1, s1)
    g2 = ROOT.RooGaussian(f"g2_{tag}", "g2", x, mu2, s2)
    g3 = ROOT.RooGaussian(f"g3_{tag}", "g3", x, mu3, s3)

    c1 = ROOT.RooRealVar(f"c1_{tag}", "c_{1}", 0.5, 0.0, 1.0)
    c2 = ROOT.RooRealVar(f"c2_{tag}", "c_{2}", 0.5, 0.0, 1.0)
    # recursive fractions: f1 = c1, f2 = (1-c1)*c2, f3 = (1-c1)*(1-c2)  -> always sum to 1
    pdf = ROOT.RooAddPdf(f"pdf_{tag}", "3 Gaussians",
                         ROOT.RooArgList(g1, g2, g3), ROOT.RooArgList(c1, c2), True)

    free = [mu1, s1, d2, d3, c1, c2] + ([] if common_mean else [mu2, mu3])
    start = {p.GetName(): p.getVal() for p in free}
    rng = np.random.default_rng(seed)

    best = None
    attempts = []
    xs_win = np.linspace(lo_fit, hi_fit, 3001)
    for k in range(n_tries):
        for p in free:
            v = start[p.GetName()]
            if k > 0:  
                span = p.getMax() - p.getMin()
                v = v + 0.15 * span * rng.uniform(-1, 1)
                v = min(max(v, p.getMin() + 1e-3 * span), p.getMax() - 1e-3 * span)
            p.setVal(v)

        res = pdf.fitTo(dh, RF.Range("fitRange"), RF.SumW2Error(True), RF.Save(True),
                        RF.PrintLevel(-1), RF.Strategy(2), RF.Minimizer("Minuit2", "migrad"))
        
        if res.status() == 0 and res.covQual() == 3:
            rank = 0
        elif res.status() in (0, 1, 2):
            rank = 1
        else:
            rank = 2
        pk_try = float(xs_win[np.argmax(scan_pdf(pdf, x, xs_win))])
        rec = dict(rank=rank, nll=res.minNll(), status=res.status(), covq=res.covQual(),
                   peak=pk_try, vals={p.GetName(): p.getVal() for p in free})
        attempts.append(rec)
        if best is None or (rec["rank"], rec["nll"]) < (best["rank"], best["nll"]):
            best = rec

    n_ok = sum(a["rank"] == 0 for a in attempts)
    n_min = sum(a["rank"] <= 1 for a in attempts)
    good = [a["peak"] for a in attempts if a["rank"] <= 1] or [a["peak"] for a in attempts]
    peak_spread = max(good) - min(good)

    for p in free:  # restore the best attempt
        p.setVal(best["vals"][p.GetName()])

    at_limit = []
    for p in free:
        span = p.getMax() - p.getMin()
        if abs(p.getVal() - p.getMin()) < 1e-2 * span or abs(p.getVal() - p.getMax()) < 1e-2 * span:
            at_limit.append(p.GetName().split("_")[0])

    # peak = maximum of the fitted pdf inside the fit window ; effSigma of the fitted pdf
    xs = np.linspace(XLO, XHI, 20001)
    vals = scan_pdf(pdf, x, xs)
    win = (xs >= lo_fit) & (xs <= hi_fit)
    peak_fit = float(xs[win][np.argmax(vals[win])])
    es_pdf = float(effsigma_from_curve(xs, vals))

    # chi2/ndf in the fit range
    fch = x.frame(RF.Range("fitRange"))
    dh.plotOn(fch, RF.DataError(ROOT.RooAbsData.SumW2), RF.Name("d"))
    pdf.plotOn(fch, RF.Range("fitRange"), RF.NormRange("fitRange"), RF.Name("m"))
    n_float = len(free)
    chi2ndf = fch.chiSquare("m", "d", n_float)
    centers = np.array([h.GetBinCenter(i) for i in range(1, nbins + 1)])
    ndf = int(((centers >= lo_fit) & (centers <= hi_fit)).sum()) - n_float
    f1 = c1.getVal()
    f2 = (1 - c1.getVal()) * c2.getVal()
    f3 = (1 - c1.getVal()) * (1 - c2.getVal())

    return dict(
        tag=tag, x=x, dh=dh, pdf=pdf, gauss=[g1, g2, g3], lo_fit=lo_fit, hi_fit=hi_fit,
        peak_fit=peak_fit, es_pdf=es_pdf, chi2ndf=chi2ndf, ndf=ndf, chi2=chi2ndf * ndf,
        status=best["status"], covq=best["covq"], rank=best["rank"],
        n_ok=n_ok, n_min=n_min, n_tries=n_tries, at_limit=at_limit, peak_spread=peak_spread,
        fracs=(f1, f2, f3),
        mus=(mu1.getVal(), mu2.getVal(), mu3.getVal()),
        sigmas=(s1.getVal(), s2.getVal(), s3.getVal()),
        _keep=[mu1, mu2, mu3, s1, d2, d3, s2, s3, c1, c2, g1, g2, g3, fch],
    )


def fit_peak_roofit(masses, weights, eff_sigma, tag, nbins=160, half_win=4.0,
                    n_tries=5, n_pass=2, common_mean=False, seed=1):
    """Histogram -> (re-centred) RooFit 3-Gaussian fit. Returns (fit dict, peak_hist, bin width)."""
    h = ROOT.TH1D(f"h_{tag}", "", nbins, XLO, XHI)
    h.Sumw2()
    m = np.ascontiguousarray(masses, dtype=np.float64)
    w = np.ones_like(m) if weights is None else np.ascontiguousarray(weights, dtype=np.float64)
    h.FillN(len(m), m, w)

    # Rescale to the effective number of entries: shape unchanged, but MINUIT sees a well-scaled NLL
    # (tiny weights, e.g. lumi/sumw ~ 1e-6, make the NLL almost flat -> bad convergence).
    neff = w.sum() ** 2 / (w ** 2).sum()
    if h.Integral() > 0:
        h.Scale(neff / h.Integral())

    peak_hist = h.GetXaxis().GetBinCenter(h.GetMaximumBin())
    center = peak_hist
    fr = None
    for ipass in range(n_pass):
        fr = do_fit(h, center, eff_sigma, f"{tag}_p{ipass}", half_win, n_tries, common_mean, seed)
        if abs(fr["peak_fit"] - center) < 0.05:
            break
        center = fr["peak_fit"]  # re-centre the +-half_win window on the fitted peak
    fr["h"] = h
    fr["neff"] = neff
    return fr, peak_hist, h.GetBinWidth(1)


# plot with all info, not all is required...clean up later
# ----------------------------------------------------------------------
def make_plot(fr, info, plot_dir, show_info=True):
    x, dh, pdf = fr["x"], fr["dh"], fr["pdf"]
    lo_fit, hi_fit, pk = fr["lo_fit"], fr["hi_fit"], fr["peak_fit"]
    es_lo, es_hi = info["win"]
    tag = fr["tag"]

    frame = x.frame(RF.Title(" "))
    dh.plotOn(frame, RF.DataError(ROOT.RooAbsData.SumW2), RF.Name("data"))
    pdf.plotOn(frame, RF.Range(XLO, XHI), RF.NormRange("fitRange"), RF.LineColor(ROOT.kGray + 1),
               RF.LineStyle(ROOT.kDotted), RF.LineWidth(2), RF.Name("ext"))
    pdf.plotOn(frame, RF.Range("fitRange"), RF.NormRange("fitRange"),
               RF.LineColor(ROOT.kBlue + 1), RF.LineWidth(3), RF.Name("model"))
    colors = [ROOT.kGreen + 2, ROOT.kMagenta + 1, ROOT.kOrange + 7]
    for i, g in enumerate(fr["gauss"]):
        pdf.plotOn(frame, RF.Components(g.GetName()), RF.Range("fitRange"), RF.NormRange("fitRange"),
                   RF.LineStyle(ROOT.kDashed), RF.LineColor(colors[i]), RF.LineWidth(2),
                   RF.Name(f"comp{i}"))

    pull = frame.pullHist("data", "ext")

    ymax = frame.GetMaximum()
    frame.SetMinimum(0.0)
    frame.SetMaximum(1.6 * ymax)
    bw = (XHI - XLO) / fr["h"].GetNbinsX()
    frame.GetYaxis().SetTitle(f"Events / ({bw:.3f} GeV)  [norm. to N_{{eff}}]")
    frame.GetYaxis().SetTitleSize(0.045)
    frame.GetYaxis().SetTitleOffset(1.1)
    frame.GetXaxis().SetLabelSize(0)
    frame.GetXaxis().SetTitleSize(0)

    c = ROOT.TCanvas(f"c_{tag}", "", 900, 800)
    pad1 = ROOT.TPad(f"p1_{tag}", "", 0, 0.27, 1, 1)
    pad2 = ROOT.TPad(f"p2_{tag}", "", 0, 0, 1, 0.27)
    pad1.SetBottomMargin(0.02)
    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.38)
    pad1.Draw()
    pad2.Draw()

    pad1.cd()
    frame.Draw()
    top = 1.6 * ymax
    box = ROOT.TBox(es_lo, 0, es_hi, 0.62 * top)
    box.SetFillColorAlpha(ROOT.kRed, 0.12)
    box.SetLineColor(0)
    box.Draw()
    l_pk = ROOT.TLine(pk, 0, pk, 0.62 * top)
    l_pk.SetLineColor(ROOT.kRed + 1)
    l_pk.SetLineWidth(2)
    l_pk.Draw()
    l_a = ROOT.TLine(lo_fit, 0, lo_fit, 0.62 * top)
    l_b = ROOT.TLine(hi_fit, 0, hi_fit, 0.62 * top)
    for ln in (l_a, l_b):
        ln.SetLineColor(ROOT.kGray + 2)
        ln.SetLineStyle(4)
        ln.SetLineWidth(2)
        ln.Draw()

    leg = ROOT.TLegend(0.11, 0.57, 0.52, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.AddEntry(frame.findObject("data"), f"Signal MC, m_{{a}} = {info['m_val']} GeV", "pe")
    leg.AddEntry(frame.findObject("model"), "Fit: sum of 3 Gaussians", "l")
    for i in range(3):
        leg.AddEntry(frame.findObject(f"comp{i}"),
                     f"G{i+1}: f={fr['fracs'][i]:.2f}, #mu={fr['mus'][i]:.2f}, "
                     f"#sigma={fr['sigmas'][i]:.2f}", "l")
    leg.AddEntry(frame.findObject("ext"), "Fit extrapolated outside range", "l")
    leg.AddEntry(l_pk, f"Fitted peak = {pk:.3f} GeV", "l")
    leg.AddEntry(box, f"68.3% window: #sigma_{{eff}} = {info['eff_sigma']:.3f} GeV", "f")
    leg.AddEntry(l_a, f"Fit range [{lo_fit:.2f}, {hi_fit:.2f}] GeV", "l")
    leg.Draw()

    if show_info:  # black-bordered info box (toggle with --no-info-box)
        pt = ROOT.TPaveText(0.555, 0.575, 0.895, 0.895, "NDC")
        pt.SetFillColor(0)
        pt.SetBorderSize(1)
        pt.SetTextAlign(12)
        pt.SetTextSize(0.026)
        pt.AddText(f"N_{{evts}} = {info['n']}  (N_{{eff}} = {info['neff']:.0f})")
        pt.AddText(f"#chi^{{2}}/ndf = {fr['chi2']:.1f}/{fr['ndf']} = {fr['chi2ndf']:.2f}")
        pt.AddText(f"Fit status = {fr['status']}, covQual = {fr['covq']}, {['OK', 'COV', 'FAIL'][fr['rank']]}")
        pt.AddText(f"Peak spread (attempts) = {fr['peak_spread']:.3f} GeV")
        pt.AddText(f"Peak (fit) = {pk:.3f} GeV")
        pt.AddText(f"Peak (hist bin) = {info['peak_hist']:.3f} GeV")
        pt.AddText(f"#sigma_{{eff}} (data) = {info['eff_sigma']:.3f} GeV")
        pt.AddText(f"#sigma_{{eff}} (fit pdf) = {fr['es_pdf']:.3f} GeV")
        pt.AddText(f"Window = [{es_lo:.2f}, {es_hi:.2f}] GeV")
        pt.Draw()

    pad2.cd()
    frame2 = x.frame(RF.Title(" "))
    frame2.addPlotable(pull, "P")
    frame2.SetMinimum(-5)
    frame2.SetMaximum(5)
    frame2.GetYaxis().SetTitle("Pull")
    frame2.GetYaxis().SetTitleSize(0.12)
    frame2.GetYaxis().SetTitleOffset(0.4)
    frame2.GetYaxis().SetLabelSize(0.10)
    frame2.GetYaxis().SetNdivisions(505)
    frame2.GetXaxis().SetTitleSize(0.14)
    frame2.GetXaxis().SetLabelSize(0.11)
    frame2.GetXaxis().SetTitleOffset(1.0)
    frame2.Draw()
    zero = ROOT.TLine(XLO, 0, XHI, 0)
    zero.SetLineStyle(2)
    zero.Draw()
    pa = ROOT.TLine(lo_fit, -5, lo_fit, 5)
    pb = ROOT.TLine(hi_fit, -5, hi_fit, 5)
    for ln in (pa, pb):
        ln.SetLineColor(ROOT.kGray + 2)
        ln.SetLineStyle(4)
        ln.Draw()

    os.makedirs(plot_dir, exist_ok=True)
    c.SaveAs(os.path.join(plot_dir, f"fit_{info['m_val']}GeV.png"))
    c.SaveAs(os.path.join(plot_dir, f"fit_{info['m_val']}GeV.pdf"))


def main():
    p = argparse.ArgumentParser(description="effSigma + RooFit 3-Gaussian peak per mass point")
    p.add_argument("-base", "--base-dir", type=str,
                   default="/eos/user/a/arnaik/Sergi_oldPresel_DataBkgSig_toTrain/signal_processedOldPresel")
    p.add_argument("-m", "--mass-col", type=str, default="mass_gggg")
    p.add_argument("-w", "--weight-col", type=str, default="weight")
    p.add_argument("--nbins", type=int, default=160)
    p.add_argument("--half-win", type=float, default=4.0, help="Fit window half-width around the peak")
    p.add_argument("--n-tries", type=int, default=5, help="Fit attempts (jittered start values)")
    p.add_argument("--n-pass", type=int, default=2, help="Re-centre the fit window on the fitted peak")
    p.add_argument("--common-mean", action="store_true", help="Share one mean between the 3 Gaussians")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--filter-mass", action="store_true", help="Restrict to 115 < mass_gggg < 135")
    p.add_argument("--apply-scaling", action="store_true", help="Apply lumi / sumw scaling")
    p.add_argument("--lumi", type=float, default=1.0)
    p.add_argument("--save-plots", action="store_true", help="Save a fit plot per mass point")
    p.add_argument("--no-info-box", action="store_true",
                   help="Hide the black-bordered info box (chi2, status, peaks, sigmas...) in the plots")
    p.add_argument("--plot-dir", type=str, default="fit_plots")
    args = p.parse_args()

    if not args.filter_mass:
        print("NOTE: --filter-mass is off: effSigma will be computed on the full mass range.")

    mass_points = list(range(15, 65, 5))
    sumw_values = [1025691, 971967, 969906, 981050, 1047791,
                   981828, 974661, 1022301, 994105, 955078] # as listed in https://gitlab.cern.ch/takumar/h4g-tools/-/blob/decoupled_from_cmssw/h4g_tools/utils/scales.py?ref_type=heads#L62
    sumw_dict = dict(zip(mass_points, sumw_values))

    results, diags = {}, {}
    for m_val in mass_points:
        folder = os.path.join(args.base_dir, f"Signal_{m_val}_GeV", "nominal")
        files = glob.glob(os.path.join(folder, "*.parquet"))
        if not files:
            print(f"Signal_{m_val}_GeV | WARNING: no parquet files in {folder}")
            continue

        cols = [args.mass_col] + ([args.weight_col] if args.weight_col else [])
        events = ak.concatenate([ak.from_parquet(f, columns=cols) for f in files])
        masses = ak.to_numpy(events[args.mass_col])
        weights = None
        if args.weight_col and args.weight_col in events.fields:
            weights = ak.to_numpy(events[args.weight_col])

        if args.filter_mass:
            mask = (masses >= XLO) & (masses <= XHI)
            masses = masses[mask]
            if weights is not None:
                weights = weights[mask]
        if len(masses) == 0:
            print(f"Signal_{m_val}_GeV | no events after filtering")
            continue

        if args.apply_scaling:
            sf = args.lumi / sumw_dict[m_val]
            weights = np.full_like(masses, sf, dtype=float) if weights is None else weights * sf

        diag = weight_diagnostics(masses, weights)
        diags[m_val] = diag
        eff_sigma, (lo, hi) = getEffSigma_unbinned(masses, weights)
        eff_sigma_unw, _ = getEffSigma_unbinned(masses, None)  # ignoring weights, cross-check

        fr, peak_hist, _ = fit_peak_roofit(
            masses, weights, eff_sigma, tag=f"{m_val}GeV", nbins=args.nbins,
            half_win=args.half_win, n_tries=args.n_tries, n_pass=args.n_pass,
            common_mean=args.common_mean, seed=args.seed)

        results[m_val] = dict(n=len(masses), neff=diag["neff"], es=eff_sigma, es_unw=eff_sigma_unw,
                              es_pdf=fr["es_pdf"], peak=fr["peak_fit"], peak_hist=peak_hist,
                              chi2ndf=fr["chi2ndf"], status=fr["status"], covq=fr["covq"],
                              quality=["OK", "COV", "FAIL"][fr["rank"]], spread=fr["peak_spread"],
                              nok=fr["n_ok"], nmin=fr["n_min"], ntries=fr["n_tries"],
                              lim=",".join(fr["at_limit"]))

        if args.save_plots:
            info = dict(m_val=m_val, n=len(masses), neff=diag["neff"], eff_sigma=eff_sigma,
                        win=(lo, hi), peak_hist=peak_hist)
            make_plot(fr, info, args.plot_dir, show_info=not args.no_info_box)

    print("\nWeight diagnostics (negative / very uneven weights can distort effSigma)") # this was added for a check, can be removed now----
    print("=" * 80)
    print(f"{'Mass Point':<12} | {'N':<8} | {'N_eff':<9} | {'w_min':<11} | {'w_max':<11} | {'#neg w':<7}")
    print("=" * 80)
    for m_val, d in diags.items():
        print(f"{m_val:<5} GeV    | {d['n']:<8} | {d['neff']:<9.1f} | {d['wmin']:<11.3e} | "
              f"{d['wmax']:<11.3e} | {d['n_neg']:<7}")
    print("=" * 80)
    #----------------------------------------------------------------------------------------------------------------------------------------

    print("\nResults")
    print("=" * 46)
    '''
    print(f"{'Mass Point':<12} | {'effSigma':<9} | {'effSig(unw)':<11} | {'effSig(pdf)':<11} | "
          f"{'Peak fit':<9} | {'Peak hist':<9} | {'spread':<7} | {'chi2/ndf':<8} | {'Fit':<5} | "
          f"{'OK/min/tries':<12} | {'at limit'}")
          '''
    print(f"{'Mass Point':<12} | {'effSigma':<9} | {'Peak fit':<9}")
    
    print("=" * 46)
    '''
    for m_val, r in results.items():
        att = f"{r['nok']}/{r['nmin']}/{r['ntries']}"
        print(f"{m_val:<5} GeV    | {r['es']:<9.4f} | {r['es_unw']:<11.4f} | {r['es_pdf']:<11.4f} | "
              f"{r['peak']:<9.4f} | {r['peak_hist']:<9.4f} | {r['spread']:<7.4f} | {r['chi2ndf']:<8.2f} | "
              f"{r['quality']:<5} | {att:<12} | {r['lim']}")
              '''
    for m_val, r in results.items():
        print(f"{m_val:<5} GeV    | {r['es']:<9.4f} | {r['peak']:<9.4f}")
    
    print("=" * 46)
    print("effSigma = sliding window (68.3%) on data (weighted); Peak fit = max of the 3-Gaussian fit. Values in GeV.")
    
    '''
    print("effSigma = sliding window on data (weighted); effSig(unw) = same, weights ignored;")
    print("effSig(pdf) = same algorithm applied to the fitted pdf (cross-check). All values in GeV.")
    print("spread = max-min of the fitted peak over the converged attempts (stability of the peak).")
    print("Fit: OK = converged, accurate covariance | COV = minimum found, covariance only approximate")
    print("     (peak is still reliable if spread is small) | FAIL = MINUIT did not converge.")
    '''


if __name__ == "__main__":
    main()
