import os, sys, re, glob
import h4g_tools.utils.loading as l

SCRIPT = 'apply_categories.py'
arg = lambda f: sys.argv[sys.argv.index(f) + 1]
year, xgb_name, nCats = arg('-y'), arg('-xgb'), arg('--nCats')
sample = 'background' if ('-bkg' in sys.argv or '--bkg' in sys.argv) else 'signal'
here = os.path.dirname(os.path.abspath(SCRIPT))
base = os.path.join(here, 'outputs', 'BDT', year, xgb_name)

#It makes use of the cats_{1}.txt obtained from running the calculate_categories.py
#cats_txt = open(os.path.join(here, 'cats', year, xgb_name, f'cats{nCats}.txt')).read()
#--------------for now..for 2 sigma-----
SUB = '2_sigma_corrected'
cats_txt = open(os.path.join(here, 'cats', year, xgb_name, SUB, f'cats{nCats}.txt')).read()
#-------------------------------------

in_cats = {int(m) for m in re.findall(r'Mass Point: (\d+) GeV', cats_txt)}
OK = [m for m in range(15, 65, 5)
      if m in in_cats and glob.glob(f'{base}/{sample}/{m}_GeV/**/*.parquet', recursive=True)]
assert OK, 'No usable mass points'
print(f'[PATCH] Mass points: {OK}')

#skip missing mass points
_orig_loadSamples = l.loadSamples
def _loadSamples(directory, *a, **k):
    try:
        return _orig_loadSamples(directory, *a, **k)
    except AssertionError:
        d = str(directory); m = re.search(r'(\d+)_GeV/?$', d)
        if not m or int(m.group(1)) in OK:
            raise
        return _orig_loadSamples(re.sub(r'\d+_GeV(/?)$', f'{OK[0]}_GeV\\1', d), *a, **k)
l.loadSamples = _loadSamples

_orig_loadPerMass = l.loadPerMass
OK_KEYS = {f'{m}_GeV' for m in OK}
l.loadPerMass = lambda *a, **k: {x: v for x, v in _orig_loadPerMass(*a, **k).items() if x in OK_KEYS}

#--------------------------------------------
src = open(os.path.join(here, SCRIPT)).read()
for old, new in [('{cwd}/../../scripts/outputs', '{cwd}/outputs'),
                 ('[f"{m}_GeV" for m in range(15,65,5)]', 'OK_LIST')]:
    assert old in src, f'pattern not found: {old}'
    src = src.replace(old, new)

#---------just for 2 sigma----------------
# ---- read cuts from / write json into the 2_sigma subfolder ----
pat = '/cats/{args.year}/{args.xgb_name}/'
assert pat in src, f'pattern not found: {pat}'
src = src.replace(pat, pat + SUB + '/')   # replaces both the read and the write
os.makedirs(os.path.join(here, 'cats', year, xgb_name, SUB), exist_ok=True)
#-------------------------------------------------------------------

g = {'__name__': '__main__', '__file__': os.path.join(here, SCRIPT), 'OK_LIST': [f'{m}_GeV' for m in OK]}
exec(compile(src, os.path.join(here, SCRIPT), 'exec'), g)
