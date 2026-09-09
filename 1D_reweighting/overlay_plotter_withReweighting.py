import os
import numpy as np
import pandas as pd
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(1111111)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetTitleBorderSize(0)

LOGY_YMIN = 1e-6

# -------------------------------------------------------
# Input parquet files
# -------------------------------------------------------
parquet_paths = {
    "skimmed": "/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/files_to_analyze/2024C_noMassCut/skimmed_files/nominal/960b74ae-9d64-11f1-9c9c-4f69b8bcbeef_Events_0-24234.parquet",
    "nc20":    "FallbackPt0001_pTgggg_1Dreweighted_v2.parquet",
}

# Medium blue (the "Event Mixing" color)
EVENT_MIX_BLUE = ROOT.TColor.GetColor("#446CCF")  # ~#5B9BD5

samples = {
    "skimmed": {"label": "Data",           "color": ROOT.kBlack},
    "nc20":    {"label": "Event Mix (nc20)", "color": EVENT_MIX_BLUE},
}

# Per-sample weight column. Data has no weight column (weight = 1 implicitly).

WEIGHT_COL_BY_SAMPLE = {
    "skimmed": None,
    "nc20": "pT_gggg_1d_weight",
    #"nc20": None,
}

# -------------------------------------------------------
# Variable config -- identical to the original plotting script,
# plus the two dR(gamma,gamma)-per-pseudoscalar variables you asked for
# (these already exist as columns: LeadPs_dR_gg / SubleadPs_dR_gg).
# -------------------------------------------------------
var_config = {
    "hist_m_gggg": {"title": "m_{4#gamma} [GeV]", "xlim": (110, 180), "logy": True, "bins": (100, 0, 500)},
    "hist_m4g": {"title": "m_{4#gamma} [GeV]", "xlim": (0, 500), "logy": True,  "bins": (100, 0, 500)},

    "hist_cos_ag": {"title": "cos(#theta^{*})", "bins":(100, -1, 1)},
    "h_cos_ag": {"title": "cos(#theta^{*})", "bins":(100, -1, 1)},

    "hist_pho1_pt": {"title": "#gamma_{1} p_{T} [GeV]",  "xlim": (0, 750), "bins": (75, 0, 750), "logy": True},
    "hist_pho2_pt": {"title": "#gamma_{2} p_{T} [GeV]", "xlim": (0, 500), "bins": (50,0,500), "logy": True},
    "hist_pho3_pt": {"title": "#gamma_{3} p_{T} [GeV]", "xlim": (0, 500), "bins": (50,0,500), "logy": True},
    "hist_pho4_pt": {"title": "#gamma_{4} p_{T} [GeV]", "xlim": (0, 500), "bins": (50,0,500), "logy": True},
    "hist_pho1_eta": {"title": "#gamma_{1} #eta", "bins": (60,-3,3)},
    "hist_pho2_eta": {"title": "#gamma_{2} #eta", "bins": (60,-3,3)},
    "hist_pho3_eta": {"title": "#gamma_{3} #eta", "bins": (60,-3,3)},
    "hist_pho4_eta": {"title": "#gamma_{4} #eta", "bins": (60,-3,3)},
    "hist_pho1_phi": {"title": "#gamma_{1} #phi", "bins": (60,-3.15,3.15)},
    "hist_pho2_phi": {"title": "#gamma_{2} #phi", "bins": (60,-3.15,3.15)},
    "hist_pho3_phi": {"title": "#gamma_{3} #phi", "bins": (60,-3.15,3.15)},
    "hist_pho4_phi": {"title": "#gamma_{4} #phi", "bins": (60,-3.15,3.15)},
    "hist_pho1_mvaID": {"title": "#gamma_{1} MVA ID", "logy": True, "bins": (50, -1.0, 1.0)},
    "hist_pho2_mvaID": {"title": "#gamma_{2} MVA ID", "logy": True, "bins": (50, -1.0, 1.0)},
    "hist_pho3_mvaID": {"title": "#gamma_{3} MVA ID", "logy": True, "bins": (50, -1.0, 1.0)},
    "hist_pho4_mvaID": {"title": "#gamma_{4} MVA ID", "logy": True, "bins": (50, -1.0, 1.0)},

    "hist_LeadPs_mass": {"title": "m_{LeadPs} [GeV]", "bins": (75, 0, 150)},
    "hist_SubleadPs_mass": {"title": "m_{SubleadPs} [GeV]", "bins": (75, 0, 150)},
    "hist_LeadPs_pt": {"title": "LeadPs p_{T} [GeV]", "logy": True, "bins": (150, 0, 300)},
    "hist_SubleadPs_pt": {"title": "SubleadPs p_{T} [GeV]", "logy": True, "bins": (150, 0, 300)},
    "hist_Ps_massDiff": {"title": "m_{LeadPs} - m_{SubleadPs} [GeV]", "bins": (75, -75, 75)},
    "hist_dM": {"title": "|m_{LeadPs} - m_{SubleadPs}| [GeV]", "bins": (75, 0, 150)},
    "hist_dR_aa": {"title": "#DeltaR(a_{1}, a_{2})", "rebin": 2, "xlim": (0, 10), "bins": (60, 0, 6)},
    "hist_dR_aa_mass_gggg": {"title": "#DeltaR(a_{1}, a_{2}) / m_{4#gamma}", "bins": (100, 0, 0.05)},

    # NEW: per-pseudoscalar photon-photon opening angle
    "hist_LeadPs_dR_gg": {"title": "#DeltaR(#gamma,#gamma)_{LeadPs}", "xlim": (0, 10), "bins": (50, 0, 5)},
    "hist_SubleadPs_dR_gg": {"title": "#DeltaR(#gamma,#gamma)_{SubleadPs}", "xlim": (0, 10), "bins": (50, 0, 5)},

    "hist_LeadPs_interMass": {"title": "(m_{Lead} - m_{hyp}) / m_{4#gamma}", "bins": (40, -0.5, 0.5)},
    "hist_SubleadPs_interMass": {"title": "(m_{Sublead} - m_{hyp}) / m_{4#gamma}", "bins": (40, -0.5, 0.5)},

    "hist_sumPt_gggg_scalar": {"title": "S_{T} (Scalar Sum p_{T}) [GeV]", "logy": True, "bins": (100, 0, 500)},
    #"hist_pT_gggg_vec": {"title": "|#vec{p}_{T}^{4#gamma}| [GeV]", "logy": True, "bins": (170, 0, 850)},
    "hist_pT_gggg_vec": {"title": "|#vec{p}_{T}^{4#gamma}| [GeV]",
                            "bins": np.array([
                                0, 10, 20, 30, 40, 50, 65, 80, 100, 125, 160, 250, 850
                            ], dtype=float),
                            "logy": True},

    # other calculated quantities
    "hist_px_gggg": {"title": "p_{x}^{4#gamma} [GeV]", "logy": True, "bins": (200, -500, 500)},
    "hist_py_gggg": {"title": "p_{y}^{4#gamma} [GeV]", "logy": True, "bins": (200, -500, 500)},
    "hist_pz_gggg": {"title": "p_{z}^{4#gamma} [GeV]", "logy": True, "bins": (200, -500, 500)},

    #"hist_MET_4g_pt": {"title": "MET 4#gamma p_{T} [GeV]", "bins": (100, -500, 0)},
    "hist_p_miss_4g_x": {"title": "p_{miss, x}^{4#gamma} [GeV]", "logy": True, "bins": (100, 0, 500)},
    "hist_p_miss_4g_y": {"title": "p_{miss, y}^{4#gamma} [GeV]", "logy": True, "bins": (100, 0, 500)},
    "hist_p_miss_4g_z": {"title": "p_{miss, z}^{4#gamma} [GeV]", "logy": True, "bins": (100, 0, 500)},
    "hist_p_miss_4g_3D": {"title": "|#vec{p}_{miss}^{3D}| [GeV]", "logy": True, "bins": (100, 0, 500)},
}

# Direct column mapping for histograms that are already in the parquet
DIRECT_COLUMN_MAP = {
    "hist_m_gggg": "mass_gggg",
    "hist_m4g": "mass_gggg",
    "hist_cos_ag": "cos_ag",
    "h_cos_ag": "cos_ag",
    "hist_pho1_pt": "pho1_pt", "hist_pho2_pt": "pho2_pt",
    "hist_pho3_pt": "pho3_pt", "hist_pho4_pt": "pho4_pt",
    "hist_pho1_eta": "pho1_eta", "hist_pho2_eta": "pho2_eta",
    "hist_pho3_eta": "pho3_eta", "hist_pho4_eta": "pho4_eta",
    "hist_pho1_phi": "pho1_phi", "hist_pho2_phi": "pho2_phi",
    "hist_pho3_phi": "pho3_phi", "hist_pho4_phi": "pho4_phi",
    "hist_pho1_mvaID": "pho1_mvaID", "hist_pho2_mvaID": "pho2_mvaID",
    "hist_pho3_mvaID": "pho3_mvaID", "hist_pho4_mvaID": "pho4_mvaID",
    "hist_LeadPs_mass": "LeadPs_mass",
    "hist_SubleadPs_mass": "SubleadPs_mass",
    "hist_LeadPs_pt": "LeadPs_pt",
    "hist_SubleadPs_pt": "SubleadPs_pt",
    "hist_Ps_massDiff": "Ps_massDiff",
    "hist_dM": "dM",
    "hist_dR_aa": "dR_aa",
    "hist_dR_aa_mass_gggg": "dR_aa_mass_gggg",
    "hist_LeadPs_dR_gg": "LeadPs_dR_gg",
    "hist_SubleadPs_dR_gg": "SubleadPs_dR_gg",
    "hist_LeadPs_interMass": "LeadPs_interMass",
    "hist_SubleadPs_interMass": "SubleadPs_interMass",
    "hist_sumPt_gggg_scalar": "sumPt_gggg",
}

# Histograms that must be computed (not present as columns)
COMPUTED_HISTS = {
    "hist_pT_gggg_vec",
    "hist_p_miss_4g_x", "hist_p_miss_4g_y", "hist_p_miss_4g_z", "hist_p_miss_4g_3D",
    "hist_px_gggg", "hist_py_gggg", "hist_pz_gggg",
}

kinematic_energy_keys = ["_pt", "_mass", "sumPt", "p_miss", "energy", "mass", "m_gggg", "m4g"] # may come handy later


def compute_derived_columns(df):
    px_total = (df["pho1_pt"] * np.cos(df["pho1_phi"]) + df["pho2_pt"] * np.cos(df["pho2_phi"]) +
                df["pho3_pt"] * np.cos(df["pho3_phi"]) + df["pho4_pt"] * np.cos(df["pho4_phi"]))
    py_total = (df["pho1_pt"] * np.sin(df["pho1_phi"]) + df["pho2_pt"] * np.sin(df["pho2_phi"]) +
                df["pho3_pt"] * np.sin(df["pho3_phi"]) + df["pho4_pt"] * np.sin(df["pho4_phi"]))
    pz_total = (df["pho1_pt"] * np.sinh(df["pho1_eta"]) + df["pho2_pt"] * np.sinh(df["pho2_eta"]) +
                df["pho3_pt"] * np.sinh(df["pho3_eta"]) + df["pho4_pt"] * np.sinh(df["pho4_eta"]))

    df["pT_gggg_vec"] = np.sqrt(px_total**2 + py_total**2)

    df["px_gggg"] = px_total
    df["py_gggg"] = py_total
    df["pz_gggg"] = pz_total

    df["p_miss_4g_x"] = -px_total
    df["p_miss_4g_y"] = -py_total
    df["p_miss_4g_z"] = -pz_total
    #df["MET_4g_pt"] = -np.sqrt(df["p_miss_4g_x"]**2 + df["p_miss_4g_y"]**2)
    df["p_miss_4g_3D"] = np.sqrt(df["p_miss_4g_x"]**2 + df["p_miss_4g_y"]**2 + df["p_miss_4g_z"]**2)
    
    return df


def get_array(df, hname):
    """Resolve the data array (and optional weights) for a given histogram name."""
    if hname in DIRECT_COLUMN_MAP:
        col = DIRECT_COLUMN_MAP[hname]
    elif hname in COMPUTED_HISTS:
        col = hname.replace("hist_", "")
    else:
        return None
    if col not in df.columns:
        return None
    return df[col].to_numpy()

# ----this function is not used and can be removed-------------------
def auto_bin_edges(data, is_energy):
    """Same automatic binning logic as the original ROOT-making script."""
    data = data[~np.isnan(data)]
    if len(data) == 0:
        return np.array([0.0, 1.0])

    if is_energy:
        max_val = np.percentile(data, 99.9)
        stop_edge = max(1.0, np.ceil(max_val) + 1.0)
        return np.arange(0.0, stop_edge, 1.0)

    q_low, q_high = np.percentile(data, [0.01, 99.99])
    if q_low == q_high:
        q_low -= 0.5
        q_high += 0.5

    unique_vals = np.unique(data)
    if len(unique_vals) <= 20:
        nbins = len(unique_vals)
        return np.linspace(unique_vals.min() - 0.5, unique_vals.max() + 0.5, nbins + 1)

    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    if iqr > 0:
        bin_w = 2 * iqr / (len(data) ** (1 / 3))
        nbins = int(np.clip((q_high - q_low) / bin_w, 50, 300))
    else:
        nbins = 200
    return np.linspace(q_low, q_high, nbins + 1)
#------------------------------------------------------------------------------


def make_th1(hname, cfg, data, weights, name):
    """Build a TH1D from a numpy array, using explicit bins if given in var_config"""
    if "bins" in cfg:
        if isinstance(cfg["bins"], np.ndarray):
            # Explicit variable-width bin edges (for pT_gggg_vec)
            edges = cfg["bins"]
            nbins = len(edges) - 1
            xmin, xmax = edges[0], edges[-1]
            h = ROOT.TH1D(name, cfg["title"], nbins, edges)
        else:
            # existing logic
            nbins, xmin, xmax = cfg["bins"]
            h = ROOT.TH1D(name, cfg["title"], nbins, xmin, xmax)
            edges = np.linspace(xmin, xmax, nbins + 1)
    else:
        is_energy = any(k in hname for k in kinematic_energy_keys)
        edges = auto_bin_edges(data, is_energy)
        nbins = len(edges) - 1
        xmin, xmax = edges[0], edges[-1]
        #h = ROOT.TH1D(name, cfg["title"], len(edges) - 1, edges)
        h = ROOT.TH1D(name, cfg["title"], nbins, xmin, xmax)

    h.Sumw2()

    counts, _ = np.histogram(data, bins=edges, weights=weights)
    if weights is not None:
        sumw2, _ = np.histogram(data, bins=edges, weights=weights**2)
    else:
        sumw2 = counts.copy()

    '''
    for i in range(len(counts)):
        h.SetBinContent(i + 1, counts[i])
        h.SetBinError(i + 1, np.sqrt(sumw2[i]))
    h.SetEntries(len(data))
    return h
    '''
    # Fill in-range bins (1 to nbins)
    for i in range(nbins):
        h.SetBinContent(i + 1, counts[i])
        h.SetBinError(i + 1, np.sqrt(sumw2[i]))

    # Compute Underflow (< xmin) and Overflow (> xmax)
    uf_mask = data < xmin
    of_mask = data > xmax

    uf_w = weights[uf_mask] if weights is not None else None
    of_w = weights[of_mask] if weights is not None else None

    uf_count = np.sum(uf_w) if weights is not None else np.sum(uf_mask)
    of_count = np.sum(of_w) if weights is not None else np.sum(of_mask)

    uf_sumw2 = np.sum(uf_w**2) if weights is not None else uf_count
    of_sumw2 = np.sum(of_w**2) if weights is not None else of_count

    # Bin 0 = Underflow, Bin nbins + 1 = Overflow
    h.SetBinContent(0, uf_count)
    h.SetBinError(0, np.sqrt(uf_sumw2))
    h.SetBinContent(nbins + 1, of_count)
    h.SetBinError(nbins + 1, np.sqrt(of_sumw2))

    h.SetEntries(len(data))
    return h


# -------------------------------------------------------
# Load parquet files, compute derived columns
# -------------------------------------------------------
dataframes = {}
for key, path in parquet_paths.items():
    print(f"Loading {key}: {path}")
    df = pd.read_parquet(path)
    df = compute_derived_columns(df)
    dataframes[key] = df
    print(f"  -> {len(df)} events")


# DEBUG: Match Data and Mixed events by event ID and save first 100 matched events
# -------------------------------------------------------

df_data = dataframes["skimmed"]
df_mix  = dataframes["nc20"]

data_debug = df_data[
    ["event_number", "pT_gggg_vec"]
].copy()

mix_debug = df_mix[
    ["event_number", "pT_gggg_vec", "pT_gggg_1d_weight"]
].copy()

data_debug = data_debug.rename(columns={
    "pT_gggg_vec": "data_pt_gggg_vec"
})

mix_debug = mix_debug.rename(columns={
    "pT_gggg_vec": "mix_pt_gggg_vec",
    "pT_gggg_1d_weight": "weight"
})

# events with a matching event ID in the mixed sample are taken
matched = pd.merge(
    data_debug,
    mix_debug,
    on="event_number",
    how="inner"
)

matched = matched.head(100)

txt_filename = "pt_gggg_vec_first100_matched.txt"

with open(txt_filename, "w") as f:

    f.write("=" * 100 + "\n")
    f.write("First 100 events: 1-D reweighting using pt_gggg_vec\n")
    f.write("=" * 100 + "\n\n")

    f.write(
        f"{'Row':>6} | "
        f"{'Event ID':>18} | "
        f"{'Data pt_gggg_vec':>20} | "
        f"{'Mix pt_gggg_vec':>20} | "
        f"{'Weight':>14}\n"
    )

    f.write("-" * 100 + "\n")

    for i, row in matched.iterrows():

        f.write(
            f"{i:6d} | "
            f"{row['event_number']:18} | "
            f"{row['data_pt_gggg_vec']:20.6f} | "
            f"{row['mix_pt_gggg_vec']:20.6f} | "
            f"{row['weight']:14.6f}\n"
        )

    f.write("=" * 100 + "\n")
    f.write("Summary")
    f.write("-" * 100 + "\n")
    f.write("Data weight: 1 (implicit)\n")
    f.write("Event Mix weight: pT_gggg_1d_weight\n")
    f.write("The weight does NOT modify pt_gggg_vec itself.\n")
    f.write("The weight modifies the contribution of the event to the histogram.\n")
    f.write("=" * 100 + "\n")

print(
    f"\nFound {len(matched)} matched events. "
    f"Saved first {len(matched)} to: {txt_filename}"
)

# DEBUG: Checking weight normalization (if integral should be 1)
mix_weights = dataframes["nc20"]["pT_gggg_1d_weight"].to_numpy()

print("\n--- Weight normalization check ---")
print(f"Number of data events   : {len(dataframes['skimmed'])}")
print(f"Number of mixed events : {len(mix_weights)}")
print(f"Sum of weights         : {np.sum(mix_weights):.6f}")
print(f"Mean weight            : {np.mean(mix_weights):.6f}")
print(f"Min weight             : {np.min(mix_weights):.6f}")
print(f"Max weight             : {np.max(mix_weights):.6f}")
print(f"Sum weights / entries  : {np.sum(mix_weights) / len(mix_weights):.6f}")
print("----------------------------------\n")

# -------------------------------------------------------
# Build all histograms in memory: hists_by_sample[key][hist_name] = TH1D
# -------------------------------------------------------
hists_by_sample = {key: {} for key in samples}
for hname, cfg in var_config.items():
    for key, df in dataframes.items():
        data = get_array(df, hname)
        if data is None:
            continue
        wcol = WEIGHT_COL_BY_SAMPLE.get(key)
        weights = df[wcol].to_numpy() if wcol else None
        h = make_th1(hname, cfg, data, weights, f"{hname}_{key}_raw")
        hists_by_sample[key][hname] = h

# -------------------------------------------------------
# Helper functions for plotting
# -------------------------------------------------------
def create_compatible_hist(source_hist, ref_hist, new_name):
    nbins = ref_hist.GetNbinsX()
    xmin = ref_hist.GetXaxis().GetXmin()
    xmax = ref_hist.GetXaxis().GetXmax()

    matched_h = ROOT.TH1D(new_name, ref_hist.GetTitle(), nbins, xmin, xmax)
    matched_h.Sumw2()

    for b in range(1, source_hist.GetNbinsX() + 1):
        center = source_hist.GetXaxis().GetBinCenter(b)
        content = source_hist.GetBinContent(b)
        error = source_hist.GetBinError(b)

        target_bin = matched_h.FindBin(center)
        if 1 <= target_bin <= nbins:
            current_val = matched_h.GetBinContent(target_bin)
            current_err = matched_h.GetBinError(target_bin)
            matched_h.SetBinContent(target_bin, current_val + content)
            matched_h.SetBinError(target_bin, (current_err**2 + error**2)**0.5)

    return matched_h


def rebin_custom(source_hist, new_name, nbins, xmin, xmax):
    new_h = ROOT.TH1D(new_name, source_hist.GetTitle(), nbins, xmin, xmax)
    new_h.Sumw2()

    for b in range(1, source_hist.GetNbinsX() + 1):
        center = source_hist.GetXaxis().GetBinCenter(b)
        content = source_hist.GetBinContent(b)
        error = source_hist.GetBinError(b)

        target_bin = new_h.FindBin(center)
        if 1 <= target_bin <= nbins:
            curr_val = new_h.GetBinContent(target_bin)
            curr_err = new_h.GetBinError(target_bin)
            new_h.SetBinContent(target_bin, curr_val + content)
            new_h.SetBinError(target_bin, (curr_err**2 + error**2)**0.5)

    return new_h


# -------------------------------------------------------
# Output setup
# -------------------------------------------------------
outfile = ROOT.TFile("Overlay_2Samples.root", "RECREATE")
os.makedirs("plots_checker_unweighted/overlay_normalized", exist_ok=True)

# -------------------------------------------------------
# Plotting loop 
# -------------------------------------------------------
for hname, cfg in var_config.items():
    hists_norm = {}
    valid = True

    for key in samples:
        h = hists_by_sample.get(key, {}).get(hname)
        if h is None or h.GetEntries() == 0:
            valid = False
            break
        h_cloned = h.Clone(f"{hname}_{key}_norm")
        if h_cloned.Integral() > 0:
            h_cloned.Scale(1.0 / h_cloned.Integral())
        hists_norm[key] = h_cloned

    # --- debug line ---
    if hname == "hist_pho1_eta" and valid:
        print("\n--- DEBUG: pho1_eta Integrals ---")
        for key, h in hists_norm.items():
            print(f"Sample [{key}] Integral (bin sum): {h.Integral():.6f} | Area (with bin widths): {h.Integral('width'):.6f}")
        print("---------------------------------\n")

    if not valid:
        continue

    print(f"Drawing: {hname}")

    canvas = ROOT.TCanvas(f"c_{hname}", hname, 1200, 1200)

    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.30, 1, 1)
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.00, 1, 0.30)

    pad1.SetBottomMargin(0.02)
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.35)

    pad1.SetRightMargin(0.28)
    pad2.SetRightMargin(0.28)

    pad1.Draw()
    pad2.Draw()

    pad1.cd()
    if cfg.get("logy", False):
        pad1.SetLogy()

    ymax = max(h.GetMaximum() for h in hists_norm.values())

    h_mix = hists_norm["nc20"]
    h_data = hists_norm["skimmed"]

    # --- Event Mixing: solid filled blue histogram ---
    h_mix.SetLineColor(samples["nc20"]["color"])
    h_mix.SetLineWidth(2)
    h_mix.SetFillColor(samples["nc20"]["color"])
    h_mix.SetFillStyle(1001)
    h_mix.SetMaximum(1.35 * ymax if not cfg.get("logy") else 10 * ymax)
    if cfg.get("logy", False):
        h_mix.SetMinimum(LOGY_YMIN)
    h_mix.GetXaxis().SetLabelSize(0)
    h_mix.GetXaxis().SetTitleSize(0)
    h_mix.GetYaxis().SetTitle("Normalized Units")
    if "xlim" in cfg:
        h_mix.GetXaxis().SetRangeUser(cfg["xlim"][0], cfg["xlim"][1])
    h_mix.Draw("HIST")

    # --- Statistical uncertainty: hatched band on top of the filled histogram ---
    h_stat = h_mix.Clone(f"{hname}_nc20_staterr")
    h_stat.SetDirectory(0)
    h_stat.SetFillColor(ROOT.kBlack)
    h_stat.SetFillStyle(3004)
    h_stat.SetLineColor(ROOT.kBlack)
    h_stat.SetMarkerSize(0)
    h_stat.Draw("E2 SAME")

    # --- Data: black points with error bars, drawn on top ---
    h_data.SetLineColor(samples["skimmed"]["color"])
    h_data.SetMarkerColor(samples["skimmed"]["color"])
    h_data.SetMarkerStyle(20)
    h_data.SetMarkerSize(0.8)
    h_data.SetLineWidth(1)
    h_data.GetXaxis().SetLabelSize(0)
    h_data.GetXaxis().SetTitleSize(0)
    if "xlim" in cfg:
        h_data.GetXaxis().SetRangeUser(cfg["xlim"][0], cfg["xlim"][1])
    h_data.Draw("E1 SAMES")

    leg = ROOT.TLegend(0.58, 0.64, 0.82, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.032)
    leg.AddEntry(h_mix, samples["nc20"]["label"], "f")
    leg.AddEntry(h_stat, "Statistical Uncertainty", "f")
    leg.AddEntry(h_data, samples["skimmed"]["label"], "lep")
    leg.Draw()

    pad1.Update()

    y_top = 0.90
    box_height = 0.22

    for key, h in hists_norm.items():
        st = h.FindObject("stats")
        if st:
            st.SetX1NDC(0.73)
            st.SetX2NDC(0.98)
            st.SetY1NDC(y_top - box_height)
            st.SetY2NDC(y_top)
            st.SetTextColor(samples[key]["color"])
            st.SetLineColor(samples[key]["color"])

            y_top -= (box_height + 0.03)

        leg = ROOT.TLegend(0.73, y_top - 0.14, 0.98, y_top)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.032)
        leg.AddEntry(h_mix, samples["nc20"]["label"], "f")
        leg.AddEntry(h_stat, "Statistical Uncertainty", "f")
        leg.AddEntry(h_data, samples["skimmed"]["label"], "lep")
        leg.Draw()

    pad1.Modified()
    pad1.Update()

    pad2.cd()

    ref_hist = hists_norm["skimmed"]
    first_ratio = True

    for key in ["nc20"]:
        raw_hist = hists_norm[key]

        if (raw_hist.GetNbinsX() != ref_hist.GetNbinsX() or
                raw_hist.GetXaxis().GetXmin() != ref_hist.GetXaxis().GetXmin() or
                raw_hist.GetXaxis().GetXmax() != ref_hist.GetXaxis().GetXmax()):

            r_num = create_compatible_hist(raw_hist, ref_hist, f"compat_{hname}_{key}")
        else:
            r_num = raw_hist.Clone(f"ratio_num_{hname}_{key}")

        r = r_num.Clone(f"ratio_{hname}_{key}")
        r.Sumw2()
        r.Divide(ref_hist)
        r.SetStats(False)
        r.SetFillStyle(0)
        r.SetLineColor(samples[key]["color"])
        r.SetMarkerColor(samples[key]["color"])
        r.SetMarkerStyle(20)
        r.SetMarkerSize(0.7)

        r.GetYaxis().SetTitle("Mix / Data")
        r.GetYaxis().CenterTitle()
        r.GetYaxis().SetNdivisions(505)
        r.GetYaxis().SetTitleSize(0.07)
        r.GetYaxis().SetTitleOffset(0.65)
        r.GetYaxis().SetLabelSize(0.06)

        r.GetXaxis().SetTitle(cfg["title"])
        r.GetXaxis().SetLabelSize(0.07)
        r.GetXaxis().SetTitleSize(0.08)
        r.GetXaxis().SetTitleOffset(1.1)

        r.SetMinimum(0.0)
        r.SetMaximum(2.0)

        if "xlim" in cfg:
            r.GetXaxis().SetRangeUser(cfg["xlim"][0], cfg["xlim"][1])

        if first_ratio:
            r.Draw("E1")
            first_ratio = False
        else:
            r.Draw("E1 SAMES")

    xmin = ref_hist.GetXaxis().GetXmin() if "xlim" not in cfg else cfg["xlim"][0]
    xmax = ref_hist.GetXaxis().GetXmax() if "xlim" not in cfg else cfg["xlim"][1]

    line = ROOT.TLine(xmin, 1.0, xmax, 1.0)
    line.SetLineStyle(2)
    line.SetLineColor(ROOT.kGray + 2)
    line.Draw()

    canvas.cd()
    canvas.Update()
    outfile.cd()
    canvas.Write()
    canvas.SaveAs(f"plots_checker_unweighted/overlay_normalized/{hname.replace('hist_', '').replace('h_', '')}.png")
    canvas.Close()

outfile.Close()

print("Plotting complete! Canvas saved to Overlay_2Samples.root and plots_checker_unweighted/overlay_normalized/")
