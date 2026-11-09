import glob, os, sys, re, runpy
import pandas as pd
import xgboost as xgb
import h4g_tools.utils.loading as l

#testing data.parquet path (it is obtained after running the 8D reweighting script)
#this path is hardcoded to /cms/cephfs/data/store/user/castells/outputs/BDT/<year>/<model-name> in input_plotter.py, therefore this wrapper is used.
MY_DATA_PARQUET = '/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/files_for_reweighting/h4g-tools/scripts/bdtIO/BDT_2024_input_weights_data.parquet'

#XGBoost for XGBoost version compatibility (crashes without this)
orig_fit = xgb.XGBClassifier.fit
def patched_fit(self, *args, **kwargs):
    eval_metric = kwargs.pop('eval_metric', None)
    early_stopping_rounds = kwargs.pop('early_stopping_rounds', None)
    params_to_set = {}
    if eval_metric is not None:
        params_to_set['eval_metric'] = eval_metric
    if early_stopping_rounds is not None:
        params_to_set['early_stopping_rounds'] = early_stopping_rounds
    if params_to_set:
        self.set_params(**params_to_set)
    return orig_fit(self, *args, **kwargs)
xgb.XGBClassifier.fit = patched_fit

#loaders: redirect ONLY the data parquet path, everything else goes to the original code
old_load = l.loadSamples
def patched_loadSamples(directory, *args, **kwargs):
    if 'BDT_2024_input_weights_data.parquet' in str(directory):
        df = pd.read_parquet(MY_DATA_PARQUET)
        branches = kwargs.get('branches', args[0] if args else None)
        if branches:
            df = df[[b for b in branches if b in df.columns]]
        print(f'[PATCH] loadSamples: data parquet {MY_DATA_PARQUET} ({len(df)} events)')
        return df
    return old_load(directory, *args, **kwargs)
l.loadSamples = patched_loadSamples

old_loadPerMass = l.loadPerMass
def patched_loadPerMass(inputs, *args, **kwargs):
    if 'BDT_2024_input_weights_data.parquet' in str(inputs):
        df = pd.read_parquet(MY_DATA_PARQUET)
        branches = kwargs.get('branches', args[0] if args else None)
        if branches:
            df = df[[b for b in branches if b in df.columns]]
        key = os.path.splitext(os.path.basename(MY_DATA_PARQUET))[0]
        print(f'[PATCH] loadPerMass: data parquet ({len(df)} events), key {key}')
        return {key: df}
    _p = str(inputs)
    if 'signal' in _p.lower() and not any(k in _p for k in ('signal', 'bkg', 'background')):
        kwargs.setdefault('era_type_arg', 'signal')
        print(f'[PATCH] loadPerMass: forcing era_type_arg=signal for {_p}')
    return old_loadPerMass(inputs, *args, **kwargs)
l.loadPerMass = patched_loadPerMass

#Scales.TransformBDT (deleted in commit 4e4d868, still called in run_bdt.py)
from h4g_tools.utils.scales import Scales as _ScalesCls
if not hasattr(_ScalesCls, 'TransformBDT'):
    def _TransformBDT(self, x, col):
        if not getattr(self, 'transform', False):
            return x
        raise NotImplementedError('--transform needs PrepareTransformation, which was removed in one of the prev commits')
    _ScalesCls.TransformBDT = _TransformBDT

#Skipping missing mass points(signal samples)
import h4g_tools.bdt.run_bdt as _rb
_orig_run_samples_BDT = _rb.run_samples_BDT
def _patched_run_samples_BDT(masses, year, path, *args, **kwargs):
    if kwargs.get('sig'):
        avail = set()
        for root, dirs, files in os.walk(str(path)):
            for name in dirs + files:
                mm = re.search(r'(\d+)\s*_?GeV', name)
                if mm:
                    avail.add(int(mm.group(1)))
        keep = [m for m in masses if m in avail]
        skipped = [m for m in masses if m not in avail]
        print(f'[PATCH] Signal masses available: {sorted(avail)}')
        print(f'[PATCH] Skipping missing masses: {skipped}')
        masses = keep
    if not kwargs.get('sig'):
        path = str(path).replace('{cwd}', os.path.dirname(os.path.abspath(_rb.__file__)))
        if 'weights_data' not in path and not os.path.exists(path):
            print(f'[PATCH] Skipping non-signal run, input not found: {path}')
            return None
    return _orig_run_samples_BDT(masses, year, path, *args, **kwargs)
_rb.run_samples_BDT = _patched_run_samples_BDT

#Plotting: only use signal masses that have saved BDT outputs
def _wrap_plot_results(fn):
    def _wrapped(masses, *args, **kwargs):
        try:
            year = args[0] if args else kwargs.get('year')
            xgb_name = args[1] if len(args) > 1 else kwargs.get('xgb_name')
            base = os.path.join(os.path.dirname(os.path.abspath(_rb.__file__)), 'outputs', 'BDT', str(year), str(xgb_name), 'signal')
            have = {int(m.group(1)) for d in os.listdir(base) for m in [re.match(r'(\d+)_GeV', d)] if m and glob.glob(os.path.join(base, d, '*', '*.parquet'))}
            keep = [m for m in masses if m in have]
            print(f'[PATCH] plot_results: masses kept {keep}, skipped {[m for m in masses if m not in have]}')
            masses = keep
        except Exception as e:
            print(f'[PATCH] plot_results mass filter not applied: {e}')
        _real = os.path.join(os.path.dirname(os.path.abspath(_rb.__file__)), 'outputs', 'BDT') + '/'
        args = tuple(_real if isinstance(a, str) and '/../../scripts/outputs' in a else a for a in args)
        print(f'[PATCH] plot_results reading from {_real}')
        _src = os.environ.get('BKG_DATA_FROM')
        if _src:
            _y = sys.argv[sys.argv.index('-y') + 1]
            _x = sys.argv[sys.argv.index('-xgb') + 1]
            _d = os.path.join(os.getcwd(), '_plot_staging', 'BDT', _y, _x)
            os.makedirs(_d, exist_ok=True)
            for _k, _n in (('signal', _x), ('background', _src), ('data', _src)):
                _l = os.path.join(_d, _k)
                if os.path.lexists(_l):
                    os.remove(_l)
                os.symlink(os.path.join(_real, _y, _n, _k), _l)
            _stage = os.path.join(os.getcwd(), '_plot_staging', 'BDT') + '/'
            args = tuple(_stage if a == _real else a for a in args)
            print(f'[PATCH] WARNING: background/data taken from {_src}, signal from {_x}. Plot-pipeline test only.')
        return fn(masses, *args, **kwargs)
    return _wrapped
import sys as _sys
for _name, _mod in list(_sys.modules.items()):
    if _name.startswith('h4g_tools') and hasattr(_mod, 'plot_results'):
        setattr(_mod, 'plot_results', _wrap_plot_results(getattr(_mod, 'plot_results')))

#skip if its hardcoded signal mass (30 GeV) has no output
_orig_generateESPlots = _rb.generateESPlots
def _patched_generateESPlots(*a, **k):
    if k.get('generic'):
        try:
            return _orig_generateESPlots(*a, **k)
        except AssertionError as e:
            print(f'[PATCH] Skipping generic plot (signal/30_GeV has no output): {e}')
            return None
    return _orig_generateESPlots(*a, **k)
_rb.generateESPlots = _patched_generateESPlots

#-------------------------------------
#adding temporarily for a check
#TEMP: train on a single signal mass point, e.g. TRAIN_SIG_MASS=30
_orig_gen_BDT = _rb.gen_BDT
def _patched_gen_BDT(*a, **k):
    m = os.environ.get('TRAIN_SIG_MASS')
    if m and k.get('eras'):
        eras = dict(k['eras'])
        eras['signal'] = [e for e in eras['signal'] if re.search(rf'_{m}_GeV(_|$)', e)]
        assert eras['signal'], f'No signal era matches mass {m}'
        print(f'[PATCH] Training with signal eras: {eras["signal"]}')
        k['eras'] = eras
    return _orig_gen_BDT(*a, **k)
_rb.gen_BDT = _patched_gen_BDT

#-------------
runpy.run_path('scan_masses_bdt.py', run_name='__main__')
