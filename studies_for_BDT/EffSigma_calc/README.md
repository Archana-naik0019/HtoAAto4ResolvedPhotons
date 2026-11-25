# Signal mass window: effSigma + 3-Gaussian peak

`2sig_win_3gauss.py` computes, for each signal mass point (15–60 GeV, steps of 5 GeV), the quantities needed to define the **signal mass window** in the H → aa → 4γ analysis:

- **effSigma**: half-width of the narrowest window containing 68.3% of the (weighted) signal events in `mass_gggg`.
- **Peak**: position of the maximum of a sum-of-3-Gaussians fit to the `mass_gggg` distribution.
- **Window**: `peak ± 2 × effSigma`.

## Process overview

```
signal parquet files (per mass point)
        │   columns: mass_gggg, weight
        ▼
optional cut 115 < mass_gggg < 135      (--filter-mass)
optional weight scaling lumi / sumw     (--apply-scaling)
        │
        ├──────────────────────────────┐
        ▼                              ▼
 effSigma                         Peak
 sliding window, 68.3%            1. TH1 of mass_gggg, 115–135 GeV, 160 bins
 (weighted)                       2. rescale histogram to N_eff
 effSigma = (hi - lo) / 2         3. first guess = highest bin
        │                         4. RooFit fit, 3 Gaussians, in [guess ± 4 GeV]
        │                         5. re-centre window on fitted peak, refit
        │                         6. peak = maximum of the fitted pdf
        ▼                              ▼
              window = peak ± 2·effSigma
                        │
                        ▼
      printed table (mass point, effSigma, peak, fit quality)
      optional fit plots with pulls    (--save-plots)
```

## Some details on the method used

**effSigma (sliding window).** The signal shape is not a pure Gaussian (it has tails), so an RMS or a single-Gaussian σ would be biased by the tails. The narrowest interval holding 68.3% of the events is a robust, fit-independent width, and it equals σ for a Gaussian. Weights are included so that the result reflects the weighted MC.

**Peak from a sum of 3 Gaussians.** A single Gaussian does not describe the core and tails of the signal, so three Gaussians are used:
- Fractions are recursive (`f1 = c1`, `f2 = (1-c1)·c2`, `f3 = (1-c1)(1-c2)`), so they sum to 1 by construction.
- Widths are ordered `σ1 < σ2 < σ3` via `σ2 = σ1 + d2`, `σ3 = σ2 + d3` with `d ≥ 0`, so the Gaussians cannot swap roles. Each width is capped so that no Gaussian becomes flat inside the fit window.
- The means are free, so the peak is not any single `μ`. It is taken as the maximum of the fitted pdf, found by a fine scan inside the fit window.

**Histogram rescaled to N_eff.** With `--apply-scaling` the weights are tiny (about lumi/sumw ~ 1e-6), which makes the NLL nearly flat and hurts MINUIT convergence. The histogram is rescaled to the effective number of entries `N_eff = (Σw)² / Σw²`. The shape is unchanged, and `SumW2Error` still gives correct uncertainties for weighted data.

**Fit window and re-centring.** The first guess from the highest bin can be off by a bin or two, so the ±4 GeV fit window is re-centred on the fitted peak and the fit is repeated until the peak moves by less than 0.05 GeV (maximum `--n-pass` passes).

**Fit quality check.** The output table includes a fit-quality flag (`OK`, `COV`, `FAIL`) so that a bad fit at any mass point does not silently shift the window.
- `OK`: converged, accurate covariance.
- `COV`: minimum found, covariance approximate (peak still reliable).
- `FAIL`: MINUIT did not converge.

## Requirements

- Python 3 with ROOT (PyROOT, including RooFit)
- `awkward`, `numpy`
- All this is ensured with the LCG stack:
```bash
#using LCG stack
source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc15-opt/setup.sh
export PYTHONPATH="$HOME/.local/lib/python3.13/site-packages/:$PYTHONPATH"
export PATH="$HOME/.local/bin/:$PATH"
```

## Usage

```bash
python3 2sig_win_3gauss.py --apply-scaling --filter-mass --save-plots
```

Expected input layout: `<base-dir>/Signal_<mass>_GeV/nominal/*.parquet`

### Options

| Option | Default | Meaning |
|---|---|---|
| `-base`, `--base-dir` | EOS signal directory | Directory containing `Signal_<m>_GeV/nominal/` |
| `-m`, `--mass-col` | `mass_gggg` | Mass column |
| `-w`, `--weight-col` | `weight` | Weight column |
| `--filter-mass` | off | Restrict to 115 < mass_gggg < 135 before computing effSigma |
| `--apply-scaling` | off | Multiply weights by `lumi / sumw` (per mass point) |
| `--lumi` | 1.0 | Luminosity used with `--apply-scaling` |
| `--nbins` | 160 | Histogram bins over 115–135 GeV |
| `--half-win` | 4.0 | Half-width of the fit window around the peak (GeV) |
| `--n-pass` | 2 | Maximum number of re-centring passes |
| `--common-mean` | off | Share one mean between the 3 Gaussians |
| `--save-plots` | off | Save a fit plot (PNG and PDF) per mass point |
| `--no-info-box` | off | Hide the info box in the plots |
| `--plot-dir` | `fit_plots` | Output directory for plots |

## Output

A table with, for each mass point, effSigma and the fitted peak (GeV), plus the fit-quality flag. The mass window used downstream is `peak ± 2 × effSigma`.

With `--save-plots`, `fit_plots/fit_<mass>GeV.{png,pdf}` show the data, the fit, the three components, the effSigma window, the fit range and the pulls.

## Notes

- `--filter-mass` matters: without it, effSigma is computed on the full mass range, including events outside 115–135 GeV.
- `--apply-scaling` does not change effSigma or the peak, since both are invariant under an overall rescaling of the weights. So it is kind of a redundant step.
- The `sumw` values used for scaling are hard-coded in the script, one per mass point.
