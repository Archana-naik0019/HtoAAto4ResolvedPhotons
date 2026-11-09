import os, sys, glob
import pandas as pd
import h4g_tools.utils.loading as l

SCRIPT = 'calculate_categories.py'
year = sys.argv[sys.argv.index('-y') + 1]
xgb_name = sys.argv[sys.argv.index('-xgb') + 1]
here = os.path.dirname(os.path.abspath(SCRIPT))
base = os.path.join(here, 'outputs', 'BDT', year, xgb_name)

#masses with parquet files in signal, background , data
def has(sample, m):
    return bool(glob.glob(f'{base}/{sample}/{m}_GeV/**/*.parquet', recursive=True))
OK = [m for m in range(15, 65, 5) if all(has(s, m) for s in ('signal', 'background', 'data'))]
OK_KEYS = {f'{m}_GeV' for m in OK}
print(f'[PATCH] Running mass points: {OK}')
print(f'[PATCH] Skipping: {[m for m in range(15, 65, 5) if m not in OK]}')
assert OK, 'No complete mass points found'

#loadSamples: skip and not crash for a missing mass folder, for testing purpose
import re
_orig_loadSamples = l.loadSamples
def _loadSamples(directory, *a, **k):
    try:
        return _orig_loadSamples(directory, *a, **k)
    except AssertionError:
        d = str(directory)
        m = re.search(r'(\d+)_GeV/?$', d)
        if not m or int(m.group(1)) in OK:
            raise
        fallback = re.sub(r'\d+_GeV(/?)$', f'{OK[0]}_GeV\\1', d)
        return _orig_loadSamples(fallback, *a, **k)
l.loadSamples = _loadSamples

_orig_loadPerMass = l.loadPerMass
def _loadPerMass(*a, **k):
    res = _orig_loadPerMass(*a, **k)
    return {key: v for key, v in res.items() if key in OK_KEYS}
l.loadPerMass = _loadPerMass



#SmoothHist: newer PyROOT returns a std::tuple proxy that can't be unpacked, so this is a jugaad..lets see
import ROOT
ROOT.gROOT.LoadMacro(os.path.join(here, '..', 'h4g_tools', 'bdt', 'smoother.cpp'))
_orig_SmoothHist = ROOT.SmoothHist
def _SmoothHist(*a, **k):
    t = _orig_SmoothHist(*a, **k)
    return (ROOT.std.get[0](t), ROOT.std.get[1](t))
ROOT.SmoothHist = _SmoothHist


#------------
src = open(os.path.join(here, SCRIPT)).read()
old = '[f"{m}_GeV" for m in range(15,65,5)]'
assert old in src, 'mass list pattern not found in script'
src = src.replace(old, 'OK_LIST')
# ---- redirect outputs into a 2_sigma subfolder ----
SUB = '2_sigma_corrected'
for pat in ['/cats/{args.xgb_name}/',                 # plots (main block)
            '/cats/{xgb_name}/',                      # plots (inside plot_ams_score)
            '/cats/{year}/{xgb_name}/',               # ams_vs_bdt.json
            '/cats/{args.year}/{args.xgb_name}/']:    # cats{n}.txt
    assert pat in src, f'output path pattern not found: {pat}'
    src = src.replace(pat, pat + SUB + '/')

for d in (f'{here}/plots/BDT/{year}/cats/{xgb_name}/{SUB}',
          f'{here}/cats/{year}/{xgb_name}/{SUB}'):
    os.makedirs(d, exist_ok=True)
#---------------------------------------
g = {'__name__': '__main__', '__file__': os.path.join(here, SCRIPT),
     'OK_LIST': [f'{m}_GeV' for m in OK]}
exec(compile(src, os.path.join(here, SCRIPT), 'exec'), g)
