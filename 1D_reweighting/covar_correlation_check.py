# Covariance/correlation is now obtained from a pairwise bivariate-Gaussian maximum-likelihood fit via Minuit2 minimizer (ROOT::Math::Minimizer).
# For each variable pair (i,j) we fit (mean_x, mean_y, sigma_x, sigma_y, rho) by minimizing the Gaussian NLL
# cov_ij = rho*sigma_x*sigma_y and
# corr_ij = rho from the fit result

import pandas as pd
import numpy as np
import ROOT
from array import array
from itertools import combinations

ROOT.gROOT.SetBatch(True)


input_file = (
    "/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/files_to_analyze/2024C_noMassCut/Evt_mix_nc20/nominal/044acbea-9e5f-11f1-ba29-98c3b8bcbeef_Events_0-484680.parquet"
)

output_file = (
    "2024C_variables_1Dhists_covCorr_Minuit_evtMixed_pTresorted.root"
)

MAX_EVENTS_PER_FIT = None


variables = [
    "pho1_mvaID",
    "pho2_mvaID",
    "pho3_mvaID",
    "pho4_mvaID",

    "LeadPs_pt",
    "SubleadPs_pt",

    "dR_aa_mass_gggg",
    "dR_aa",

    "LeadPs_interMass",
    "SubleadPs_interMass",

    "LeadPs_mass",
    "SubleadPs_mass",

    "Ps_massDiff",
    "cos_ag",

    "pho1_pt",
    "pho2_pt",
    "pho3_pt",
    "pho4_pt",

    "mass_gggg",
]


print("Reading parquet file...")
df = pd.read_parquet(input_file)

print(f"Number of events: {len(df)}")

# Calculate pT_gggg_vec / MET -----------------------------------------

px_total = (
    df["pho1_pt"] * np.cos(df["pho1_phi"])
    + df["pho2_pt"] * np.cos(df["pho2_phi"])
    + df["pho3_pt"] * np.cos(df["pho3_phi"])
    + df["pho4_pt"] * np.cos(df["pho4_phi"])
)

py_total = (
    df["pho1_pt"] * np.sin(df["pho1_phi"])
    + df["pho2_pt"] * np.sin(df["pho2_phi"])
    + df["pho3_pt"] * np.sin(df["pho3_phi"])
    + df["pho4_pt"] * np.sin(df["pho4_phi"])
)

pz_total = (
    df["pho1_pt"] * np.sinh(df["pho1_eta"])
    + df["pho2_pt"] * np.sinh(df["pho2_eta"])
    + df["pho3_pt"] * np.sinh(df["pho3_eta"])
    + df["pho4_pt"] * np.sinh(df["pho4_eta"])
)

df["pT_gggg_vec"] = np.sqrt(px_total**2 + py_total**2)

print(
    f"pT_gggg_vec range: "
    f"{df['pT_gggg_vec'].min():.3f} - "
    f"{df['pT_gggg_vec'].max():.3f} GeV"
)

met_px = -px_total
met_py = -py_total

df["MET"] = np.sqrt(met_px**2 + met_py**2)
df["MET_phi"] = np.arctan2(met_py, met_px)

print(
    f"MET range: "
    f"{df['MET'].min():.3f} - "
    f"{df['MET'].max():.3f} GeV"
)


def delta_phi(phi1, phi2):
    return np.abs((phi1 - phi2 + np.pi) % (2.0 * np.pi) - np.pi)


for i in range(1, 5):
    df[f"deltaPhi_pho{i}_MET"] = delta_phi(
        df[f"pho{i}_phi"],
        df["MET_phi"]
    )

df["min_deltaPhi_pho_MET"] = df[
    [
        "deltaPhi_pho1_MET",
        "deltaPhi_pho2_MET",
        "deltaPhi_pho3_MET",
        "deltaPhi_pho4_MET",
    ]
].min(axis=1)

deltaPhi_variables = [
    "deltaPhi_pho1_MET",
    "deltaPhi_pho2_MET",
    "deltaPhi_pho3_MET",
    "deltaPhi_pho4_MET",
    "min_deltaPhi_pho_MET",
]

print("\nDeltaPhi ranges:")
for var in deltaPhi_variables:
    print(f"  {var}: {df[var].min():.4f} - {df[var].max():.4f}")


# Photon pair angles & Photon to 4-photon system angles --------

ANGLE_PAIRS = list(combinations(range(1, 5), 2))  # (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)


def compute_all_3d_angles(df, pairs):
    
    pho_pt = [df[f"pho{i}_pt"].to_numpy() for i in range(1, 5)]
    pho_eta = [df[f"pho{i}_eta"].to_numpy() for i in range(1, 5)]
    pho_phi = [df[f"pho{i}_phi"].to_numpy() for i in range(1, 5)]

    n = len(df)
    angle_cols = {f"angle_pho{i}_pho{j}": np.zeros(n, dtype=float) for i, j in pairs}
    for i in range(1, 5):
        angle_cols[f"angle_pho{i}_sumPho"] = np.zeros(n, dtype=float)
    angle_cols["min_angle_pho_sumPho"] = np.zeros(n, dtype=float)

    for row in range(n):
        pho_vecs = []
        total_vec = ROOT.TVector3(0.0, 0.0, 0.0)

        for i in range(4):
            pt = pho_pt[i][row]
            eta = pho_eta[i][row]
            phi = pho_phi[i][row]

            px = pt * np.cos(phi)
            py = pt * np.sin(phi)
            pz = pt * np.sinh(eta)

            v = ROOT.TVector3(px, py, pz)
            pho_vecs.append(v)
            total_vec += v

        # Pairwise photon angles
        for (i, j) in pairs:
            angle_cols[f"angle_pho{i}_pho{j}"][row] = pho_vecs[i - 1].Angle(pho_vecs[j - 1])

        # Photon-to-sumPho angles
        sum_angles = [v.Angle(total_vec) for v in pho_vecs]
        for i in range(4):
            angle_cols[f"angle_pho{i+1}_sumPho"][row] = sum_angles[i]
        angle_cols["min_angle_pho_sumPho"][row] = min(sum_angles)

    for name, vals in angle_cols.items():
        df[name] = vals

    return df


print()
print(f"Computing 3D opening angles for {len(df)} events "
      f"-- this can take a little while for large samples...")
df = compute_all_3d_angles(df, ANGLE_PAIRS)

angle_variables = [f"angle_pho{i}_pho{j}" for i, j in ANGLE_PAIRS]

sumPho_angle_variables = [
    "angle_pho1_sumPho",
    "angle_pho2_sumPho",
    "angle_pho3_sumPho",
    "angle_pho4_sumPho",
    "min_angle_pho_sumPho",
]

print("\nPairwise Angle ranges:")
for var in angle_variables:
    print(f"  {var}: {df[var].min():.4f} - {df[var].max():.4f}")

print("\nPhoton-System Angle ranges:")
for var in sumPho_angle_variables:
    print(f"  {var}: {df[var].min():.4f} - {df[var].max():.4f}")


ALL_VARS = (
    variables
    + ["pT_gggg_vec", "MET", "MET_phi"]
    + angle_variables
    + sumPho_angle_variables
    + deltaPhi_variables
)


outfile = ROOT.TFile(output_file, "RECREATE")

# Hist ranges---------------------------------------------------

ranges = {
    "pho1_mvaID": (-1.0, 1.0),
    "pho2_mvaID": (-1.0, 1.0),
    "pho3_mvaID": (-1.0, 1.0),
    "pho4_mvaID": (-1.0, 1.0),

    "LeadPs_pt": (0.0, 500.0),
    "SubleadPs_pt": (0.0, 500.0),

    "pho1_pt": (0.0, 800.0),
    "pho2_pt": (0.0, 400.0),
    "pho3_pt": (0.0, 300.0),
    "pho4_pt": (0.0, 100.0),

    "dR_aa_mass_gggg": (0.0, 0.1),
    "dR_aa": (0.0, 10.0),

    "LeadPs_interMass": (0.0, 1.0),
    "SubleadPs_interMass": (0.0, 1.0),

    "LeadPs_mass": (0.0, 100.0),
    "SubleadPs_mass": (0.0, 100.0),

    "Ps_massDiff": (-50.0, 50.0),

    "cos_ag": (-1.0, 1.0),

    "mass_gggg": (110, 180),

    "pT_gggg_vec": (0.0, 700.0),
    "MET": (0.0, 700.0),
    "MET_phi": (-np.pi, np.pi),
}

for var in angle_variables + sumPho_angle_variables + deltaPhi_variables:
    ranges[var] = (0.0, np.pi)

#---------------------------------------
VAR_LABELS = {
    "pho1_mvaID": "#gamma_{1} MVA ID",
    "pho2_mvaID": "#gamma_{2} MVA ID",
    "pho3_mvaID": "#gamma_{3} MVA ID",
    "pho4_mvaID": "#gamma_{4} MVA ID",

    "LeadPs_pt": "p_{T}^{LeadPs} [GeV]",
    "SubleadPs_pt": "p_{T}^{SubleadPs} [GeV]",

    "dR_aa_mass_gggg": "#DeltaR(a_{1},a_{2}) / m_{4#gamma}",
    "dR_aa": "#DeltaR(a_{1},a_{2})",

    "LeadPs_interMass": "(m_{Lead}-m_{hyp}) / m_{4#gamma}",
    "SubleadPs_interMass": "(m_{Sublead}-m_{hyp}) / m_{4#gamma}",

    "LeadPs_mass": "m_{LeadPs} [GeV]",
    "SubleadPs_mass": "m_{SubleadPs} [GeV]",

    "Ps_massDiff": "m_{LeadPs} - m_{SubleadPs} [GeV]",

    "cos_ag": "cos(#theta^{*})",

    "pho1_pt": "#gamma_{1} p_{T} [GeV]",
    "pho2_pt": "#gamma_{2} p_{T} [GeV]",
    "pho3_pt": "#gamma_{3} p_{T} [GeV]",
    "pho4_pt": "#gamma_{4} p_{T} [GeV]",

    "mass_gggg": "m_{4#gamma} [GeV]",

    "pT_gggg_vec": "|#vec{p}_{T}^{4#gamma}| [GeV]",
    "MET": "E_{T}^{miss} [GeV]",
    "MET_phi": "#phi(E_{T}^{miss})",

    "deltaPhi_pho1_MET": "#Delta#phi(#gamma_{1},E_{T}^{miss})",
    "deltaPhi_pho2_MET": "#Delta#phi(#gamma_{2},E_{T}^{miss})",
    "deltaPhi_pho3_MET": "#Delta#phi(#gamma_{3},E_{T}^{miss})",
    "deltaPhi_pho4_MET": "#Delta#phi(#gamma_{4},E_{T}^{miss})",
    "min_deltaPhi_pho_MET": "#Delta#phi_{min}(#gamma,E_{T}^{miss})",

    "angle_pho1_sumPho": "#theta(#gamma_{1},#vec{p}_{4#gamma}) [rad]",
    "angle_pho2_sumPho": "#theta(#gamma_{2},#vec{p}_{4#gamma}) [rad]",
    "angle_pho3_sumPho": "#theta(#gamma_{3},#vec{p}_{4#gamma}) [rad]",
    "angle_pho4_sumPho": "#theta(#gamma_{4},#vec{p}_{4#gamma}) [rad]",
    "min_angle_pho_sumPho": "min #theta(#gamma_{i},#vec{p}_{4#gamma}) [rad]",
}

# photon-pair angles, generated programmatically for all 6 pairs
for _i, _j in ANGLE_PAIRS:
    VAR_LABELS[f"angle_pho{_i}_pho{_j}"] = f"#theta(#gamma_{{{_i}}},#gamma_{{{_j}}}) [rad]"


def label_for(var):
    """Symbol label if we have one, otherwise fall back to the raw column name."""
    return VAR_LABELS.get(var, var)


nBins1D = 100


# 1D histograms for vars

print("\nCreating 1D histograms...")

for var in ALL_VARS:

    if var not in df.columns:
        print(f"WARNING: {var} not found -- skipping")
        continue

    values = (
        df[var]
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
        .to_numpy()
    )

    if var in ranges:
        vMin, vMax = ranges[var]
    else:
        vMin = np.nanmin(values)
        vMax = np.nanmax(values)

        if vMin == vMax:
            vMin -= 1
            vMax += 1

    hist1D = ROOT.TH1D(
        f"h_{var}",
        f"{var};{var};Events",
        nBins1D,
        vMin,
        vMax
    )

    for value in values:
        hist1D.Fill(value)

    hist1D.Write()

    print(
        f"  {var}: "
        f"Entries = {hist1D.GetEntries():.0f}, "
        f"Range = [{vMin}, {vMax}]"
    )


# Using minuit to get cov and correlation:

print("\nComputing covariance / correlation matrices via Minuit2 bivariate-Gaussian fits...")

covVars = [v for v in ALL_VARS if v in df.columns]

covDf = df[covVars].replace([np.inf, -np.inf], np.nan).dropna()
n_full = len(covDf)

print(f"  Events available after cleaning: {n_full}")

if MAX_EVENTS_PER_FIT is not None and n_full > MAX_EVENTS_PER_FIT:
    fit_df = covDf.sample(n=MAX_EVENTS_PER_FIT, random_state=42)
    print(f"  Subsampling to {MAX_EVENTS_PER_FIT} events per fit (MAX_EVENTS_PER_FIT set)")
else:
    fit_df = covDf

n_fit = len(fit_df)

nVars = len(covVars)

col_arrays = {v: fit_df[v].to_numpy(dtype=np.float64) for v in covVars}

means = {}
sigmas = {}
for v in covVars:
    buf = array('d', col_arrays[v])
    means[v] = ROOT.TMath.Mean(n_fit, buf)
    sigmas[v] = ROOT.TMath.StdDev(n_fit, buf)


def make_nll(x_arr, y_arr, n):
    def nll(par):
        mu_x, mu_y, sx, sy, rho = par[0], par[1], par[2], par[3], par[4]
        sx = max(sx, 1e-9)
        sy = max(sy, 1e-9)
        rho = min(max(rho, -0.999), 0.999)

        dx = (x_arr - mu_x) / sx
        dy = (y_arr - mu_y) / sy
        one_minus_rho2 = 1.0 - rho * rho

        z = dx * dx - 2.0 * rho * dx * dy + dy * dy
        quad_term = np.sum(z) / (2.0 * one_minus_rho2)
        log_norm = n * np.log(2.0 * np.pi * sx * sy * np.sqrt(one_minus_rho2))

        return float(quad_term + log_norm)
    return nll


def fit_bivariate_gaussian(x_arr, y_arr, mu_x0, mu_y0, sx0, sy0, rho0):
    n = len(x_arr)
    nll = make_nll(x_arr, y_arr, n)
    functor = ROOT.Math.Functor(nll, 5)

    minimizer = ROOT.Math.Factory.CreateMinimizer("Minuit2", "Migrad")
    minimizer.SetPrintLevel(0)
    minimizer.SetMaxFunctionCalls(5000)
    minimizer.SetMaxIterations(5000)
    minimizer.SetTolerance(0.01)
    minimizer.SetFunction(functor)

    step_x = max(abs(mu_x0), 1.0) * 0.01
    step_y = max(abs(mu_y0), 1.0) * 0.01

    minimizer.SetVariable(0, "mu_x", mu_x0, step_x)
    minimizer.SetVariable(1, "mu_y", mu_y0, step_y)
    minimizer.SetLimitedVariable(2, "sigma_x", max(sx0, 1e-6), 0.05 * max(sx0, 1e-6), 1e-9, 1e9)
    minimizer.SetLimitedVariable(3, "sigma_y", max(sy0, 1e-6), 0.05 * max(sy0, 1e-6), 1e-9, 1e9)
    minimizer.SetLimitedVariable(4, "rho", rho0, 0.05, -0.999, 0.999)

    ok = minimizer.Minimize()
    result = minimizer.X()
    sx_fit, sy_fit, rho_fit = result[2], result[3], result[4]

    return bool(ok), sx_fit, sy_fit, rho_fit


covMatrix = np.zeros((nVars, nVars))
corrMatrix = np.zeros((nVars, nVars))

n_failed = 0

print(f"  Fitting {nVars * (nVars - 1) // 2} off-diagonal pairs with Minuit2 "
      f"({n_fit} events/fit)...")

for i in range(nVars):
    vi = covVars[i]

    covMatrix[i, i] = sigmas[vi] ** 2
    corrMatrix[i, i] = 1.0

    for j in range(i + 1, nVars):
        vj = covVars[j]

        x_arr = col_arrays[vi]
        y_arr = col_arrays[vj]

        cx = x_arr - means[vi]
        cy = y_arr - means[vj]
        sij0 = sigmas[vi] * sigmas[vj]
        rho0 = (np.sum(cx * cy) / (n_fit - 1)) / sij0 if sij0 != 0 else 0.0
        rho0 = min(max(rho0, -0.9), 0.9)

        ok, sx_fit, sy_fit, rho_fit = fit_bivariate_gaussian(
            x_arr, y_arr, means[vi], means[vj], sigmas[vi], sigmas[vj], rho0
        )

        if not ok:
            n_failed += 1
            print(f"    Minuit fit did not converge cleanly for ({vi}, {vj}) "
                  f"-- using its result anyway")

        cov_ij = rho_fit * sx_fit * sy_fit

        covMatrix[i, j] = cov_ij
        covMatrix[j, i] = cov_ij
        corrMatrix[i, j] = rho_fit
        corrMatrix[j, i] = rho_fit

print(f"  Done. {n_failed} / {nVars * (nVars - 1) // 2} pairs flagged as non-converged.")

hCov = ROOT.TH2D(
    "h_covMatrix",
    "Covariance Matrix (Minuit2 bivariate-Gaussian fits)",
    nVars, 0, nVars,
    nVars, 0, nVars
)

hCorr = ROOT.TH2D(
    "h_corrMatrix",
    "Correlation Matrix (Minuit2 bivariate-Gaussian fits)",
    nVars, 0, nVars,
    nVars, 0, nVars
)

for i, vname in enumerate(covVars):
    vlabel = label_for(vname)
    hCov.GetXaxis().SetBinLabel(i + 1, vlabel)
    hCov.GetYaxis().SetBinLabel(i + 1, vlabel)
    hCorr.GetXaxis().SetBinLabel(i + 1, vlabel)
    hCorr.GetYaxis().SetBinLabel(i + 1, vlabel)

for i in range(nVars):
    for j in range(nVars):
        hCov.SetBinContent(i + 1, j + 1, covMatrix[i, j])
        hCorr.SetBinContent(i + 1, j + 1, corrMatrix[i, j])

hCov.Write()
hCorr.Write()

print(f"  Wrote covariance matrix  -> h_covMatrix ({nVars}x{nVars})")
print(f"  Wrote correlation matrix -> h_corrMatrix ({nVars}x{nVars})")

def draw_matrix_png(hist, outname, zmin=None, zmax=None):
    c = ROOT.TCanvas(f"c_{outname}", hist.GetTitle(), 1600, 1400)
    c.SetLeftMargin(0.18)
    c.SetBottomMargin(0.20)
    c.SetRightMargin(0.14)
    c.SetTopMargin(0.08)

    hist.SetStats(0)
    hist.GetXaxis().SetLabelSize(0.018)
    hist.GetYaxis().SetLabelSize(0.018)
    hist.GetXaxis().LabelsOption("v")
    if zmin is not None and zmax is not None:
        hist.GetZaxis().SetRangeUser(zmin, zmax)

    hist.Draw("COLZ")
    c.SaveAs(outname)
    c.Close()
    print(f"  -> wrote {outname}")


draw_matrix_png(hCov, "covariance_matrix.png")
draw_matrix_png(hCorr, "correlation_matrix.png", zmin=-1.0, zmax=1.0)

#============================================

outfile.Close()

print("\n==========================================")
print("Finished")
print(f"Output ROOT file:")
print(output_file)
