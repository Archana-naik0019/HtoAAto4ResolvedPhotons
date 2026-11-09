import pandas as pd
import numpy as np
import xgboost
from typing import Tuple, Optional, List, Dict
import time
import os
import json
#import wandb
#from wandb.xgboost import WandbCallback
from hyperopt import fmin, tpe, hp, STATUS_OK
from sklearn.metrics import roc_auc_score, log_loss, f1_score
import h4g_tools.utils.loading as loading
import h4g_tools.utils.plotting as plotting
import h4g_tools.bdt.input_plotter as input_plotter

import sys
from sklearn.metrics import roc_curve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def define_variables(
    samplePath: str,
    reduce_df: Optional[bool] = True,
    extra: Optional[List[str]] = None,
    mass: bool = False,
    perMass: bool = False,
    eras: list = None,
    variation: str = None,
    applyReduce: bool = False,
    year: str = None
) -> pd.DataFrame:
    """
    Keep/define variables to use as inputs/training.
    """

    # Keep only training variables
    training = [
        "pho1_mvaID",
        "pho2_mvaID",
        "pho3_mvaID",
        "pho4_mvaID",
        "LeadPs_pt",
        "SubleadPs_pt",
        "LeadPs_mass",  # These are removed by reduce dataframe before training the BDT!
        "SubleadPs_mass",
        "dR_aa_mass_gggg",
        "LeadPs_interMass",
        "SubleadPs_interMass",
        "Ps_massDiff",
        "cos_ag",
        #"pT1_ma1",
        #"pT2_ma1",
        #"pT1_ma2",
        #"pT2_ma2",
    ]

    if mass:
        training.append("mass_gggg")

    if extra is not None:
        training.extend(extra)

    if reduce_df:
        dataset = loading.loadSamples(samplePath, branches=training, eras = eras, variation=variation) if not perMass else loading.loadPerMass(samplePath, branches=training, eras=eras, variation=variation)
        print(f"Kept {len(dataset.columns if not perMass else dataset[list(dataset.keys())[0]].columns)} columns.")
    else:
        dataset = loading.loadSamples(samplePath, eras = eras, variation=variation) if not perMass else loading.loadPerMass(samplePath, eras=eras, variation=variation)
        print(f"Kept {len(dataset.columns if not perMass else dataset[list(dataset.keys())[0]].columns)} columns.")

    new_dataset = {}
    if applyReduce:
        for mass in dataset.keys():
            new_dataset[int(mass[mass.find("GeV")-3 : mass.find("GeV")-1])] = reduce_dataframe(dataset[mass], extra=extra)
    else:
        new_dataset = dataset

    return new_dataset

def redefine_dataframe(
    dataset: pd.DataFrame,
    hypMass: float,
    rand: bool = False
) -> pd.DataFrame:
    """
    Corrects columns with m_Hyp for the desired new value of m_Hyp.
    """
    
    # Add custom variables for training
    dataset_updated = dataset.copy()
    if rand:
        if len(hypMass) > len(dataset):
            hypMass = hypMass[:len(dataset)]
        dataset_updated.loc[:, "LeadPs_interMass"] = (dataset.LeadPs_mass - hypMass) / dataset.mass_gggg
        dataset_updated.loc[:, "SubleadPs_interMass"] = (dataset.SubleadPs_mass - hypMass) / dataset.mass_gggg
    else:
        dataset_updated.loc[:, "LeadPs_interMass"] = (dataset.LeadPs_mass - np.full(dataset.LeadPs_mass.shape, hypMass)) / dataset.mass_gggg
        dataset_updated.loc[:, "SubleadPs_interMass"] = (dataset.SubleadPs_mass - np.full(dataset.SubleadPs_mass.shape, hypMass)) / dataset.mass_gggg

    print("Changed 2 columns for new m_Hyp.")
    
    return dataset_updated

def reduce_dataframe(
    df: pd.DataFrame,
    extra: List = None
) -> pd.DataFrame:
    """
    Reduces size of DataFrame to the columns needed for BDT training.
    """

    # Keep only training variables
    training = [
        "pho1_mvaID",
        "pho2_mvaID",
        "pho3_mvaID",
        "pho4_mvaID",
        "LeadPs_pt",
        "SubleadPs_pt",
        "dR_aa_mass_gggg",
        "LeadPs_interMass",
        "SubleadPs_interMass",
        "Ps_massDiff",
        "cos_ag",
        #"pT1_ma1",
        #"pT2_ma1",
        #"pT1_ma2",
        #"pT2_ma2",
    ]

    if extra is not None:
        training.extend(extra)

    return df[training]

TEST_MOD = 5 # added for the 20-80 split
def split_dataset(
    dataset: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split dataset into training and testing. Use even and odd rows.
    """

    #test = dataset.iloc[lambda x: x.index % 2 == 0]  # Even rows only
    #train = dataset.iloc[lambda x: x.index % 2 == 1]  # Odd rows only # changing this to check how it goes

    is_test = (dataset.index % TEST_MOD == 0)
    test = dataset[is_test]
    train = dataset[~is_test]
    
    #train = dataset.head(int(len(dataset)*0.7))  # 70% used for training
    #test = dataset.head(int(len(dataset)*0.3))  # 30% used for testing

    return (train, test)


def prep_samples(
    signal: pd.DataFrame,
    bkg_weights: pd.DataFrame,
    year: str,
    Scales: str = None
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Prepare samples for BDT training.
    """

    # Add categories to signal/background
    bkg_weights.insert(len(bkg_weights.columns), "category", [0.0] * len(bkg_weights))
    signal.insert(len(signal.columns), "category", [1.0] * len(signal))

    # Drop events with large weights to avoid sub-optimal training (reconfirmed this on 22 Apr 2026)
    drop_rate = len(signal)
    signal = signal[signal.weight < 10.0]
    drop_rate = (1 - len(signal) / drop_rate) * 100.0
    print(f"Dropped {drop_rate:.3}% of signal events due to large event weights.")

    # ---- Per-mass train/test counts (same odd/even rule as split_dataset) ----************************************
    mpt = ((signal["LeadPs_mass"] - signal["LeadPs_interMass"] * signal["mass_gggg"]) / 5).round() * 5
    #tr_m = mpt[signal.index % 2 == 1].value_counts().sort_index()
    #te_m = mpt[signal.index % 2 == 0].value_counts().sort_index()
    tr_m = mpt[signal.index % TEST_MOD != 0].value_counts().sort_index()
    te_m = mpt[signal.index % TEST_MOD == 0].value_counts().sort_index()
    print(pd.DataFrame({"train": tr_m, "test": te_m, "train/test": tr_m / te_m}))
    #------------------------**************************************************************************************

    # NOTE: weights are created from genWeight and then modified for systematics and/or corrections in processor, so this is logically equivalent to scaling based on genWeight

    # Split, set weights, and shuffle datasets
    assert "Ndimreweight" in bkg_weights.columns
    # NOTE on bkg_weights: it contains events AND weights so order is consistent.
    bkg_train, bkg_test = split_dataset(bkg_weights)

    sig_train, sig_test = split_dataset(signal)

    # Rename Ndimreweight to weight so training weight column naming is consistent
    bkg_train.rename(columns={"Ndimreweight": "weight"}, inplace=True)
    bkg_test.rename(columns={"Ndimreweight": "weight"}, inplace=True)
    # Try normalizing bkg samples to integral of sig_train at Nancy's request
    print("~~~ NORMALIZING BKG TO SIG ~~~")
    # NOTE: This is what gets my model consistent with what was trained before (even with moving sample_weight around!)
    # The weights must be balanced in some way. Doesn't matter which way...
    bkg_train["weight"] = sum(sig_train["weight"]) / sum(bkg_train["weight"]) * bkg_train["weight"]
    bkg_test["weight"] = sum(sig_train["weight"]) / sum(bkg_train["weight"]) * bkg_test["weight"]

    # Normalization of bkg to data is done in the loading of the Ndimreweights. Useful so everything is scaled to the correct statistcs.

    dtrain = pd.concat([sig_train.reset_index(drop=True), bkg_train.reset_index(drop=True)], copy=True, ignore_index=True, axis=0)
    dtest = pd.concat([sig_test.reset_index(drop=True), bkg_test.reset_index(drop=True)], copy=True, ignore_index=True, axis=0)
    dtrain = dtrain.sample(frac=1.0).reset_index(drop=True)  # Shuffle
    dtest = dtest.sample(frac=1.0).reset_index(drop=True)  # Shuffle

    return sig_train, bkg_train, reduce_dataframe(dtrain, ["category", "mass_gggg", "weight"]), reduce_dataframe(dtest, ["category", "mass_gggg", "weight"])


def train_ES_bdt(
    signal: pd.DataFrame,
    bkg_weights: pd.DataFrame,
    gen: str,
    year: str,
    modelPath: Optional[str] = None,
    save: Optional[bool] = False,
    debug: Optional[bool] = False,
    train_no_params: Optional[bool] = False,
    Scales: str = None,
    mass_gggg: dict = None,
) -> xgboost.XGBClassifier:
    """
    Setup for BDT. Calculates and add scores to events list.
    """

    start_time = time.time()

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))
    
    # Set model name for file saving
    xgb_name = modelPath.split("/")[-1].replace(".xgb","")

    # ---- DEBUG A: raw inputs, before any processing ----**********************
    DEBUG_EXIT = False   # set False once the counts look right
    m_gggg = np.asarray(mass_gggg["signal"])
    print("\n[DBG-A] signal rows:", len(signal), "| mass_gggg(sig) len:", len(m_gggg))
    print("[DBG-A] bkg rows:", len(bkg_weights), "| mass_gggg(bkg) len:", len(mass_gggg["background"]))
    print("[DBG-A] signal index unique:", signal.index.is_unique,
          "| clean 0..N-1:", signal.index.equals(pd.RangeIndex(len(signal))))
    print("[DBG-A] signal columns:", list(signal.columns))
    print("[DBG-A] variation-like columns:",
          [c for c in signal.columns if "var" in c.lower() or "syst" in c.lower()])

    mpt_raw = pd.Series(((signal["LeadPs_mass"].to_numpy()
                          - signal["LeadPs_interMass"].to_numpy() * m_gggg[:len(signal)]) / 5).round() * 5,
                        index=signal.index)
    raw = pd.DataFrame({"n_rows": mpt_raw.value_counts().sort_index()})
    raw["n_weight<10"] = signal["weight"].lt(10.0).groupby(mpt_raw).sum()
    print(raw)

    # duplicate / replication checks
    key = [c for c in ["pho1_mvaID", "pho2_mvaID", "pho3_mvaID", "pho4_mvaID",
                       "LeadPs_pt", "SubleadPs_pt"] if c in signal.columns]
    tmp = signal[key].copy()
    tmp["mpt"] = mpt_raw.to_numpy()
    tmp["w"] = signal["weight"].to_numpy()
    tmp["inter"] = signal["LeadPs_interMass"].to_numpy()
    g = tmp.groupby(key + ["mpt"])
    print("[DBG-A] unique (event, mass point):", g.ngroups, "of", len(tmp))
    print("[DBG-A] copies per event:\n", g.size().value_counts().head())
    print("[DBG-A] distinct weights per event:\n", g["w"].nunique().value_counts().head())
    print("[DBG-A] distinct interMass per event:\n", g["inter"].nunique().value_counts().head())
    if DEBUG_EXIT:
        sys.exit()
    # ---- END DEBUG A ----***************************************************

    # Prepare samples for training
    signal = pd.concat([signal, pd.Series(mass_gggg["signal"], name="mass_gggg")], axis=1)
    bkg_weights = pd.concat([bkg_weights, pd.Series(mass_gggg["background"], name="mass_gggg")], axis=1)
    sig_train, bkg_train, dtrain, dtest = prep_samples(signal, bkg_weights, year, Scales=Scales)

    # ---- DEBUG B ----****************************************
    print("\n[DBG-B] sig_train:", len(sig_train), "| bkg_train:", len(bkg_train),
          "| dtrain:", len(dtrain), "| dtest:", len(dtest))
    print("[DBG-B] NaNs in mass_gggg (signal):", signal["mass_gggg"].isna().sum())
    print("[DBG-B] sig_train index unique:", sig_train.index.is_unique)
    #sys.exit() # adding temporarily
    # ---- END DEBUG B ----***********************************

    # Plot samples here to confirm what goes into BDT exactly
    # To confirm with weights what is being used for the BDT models!
    input_plotter.plotInputs(sig_train, bkg_train, year)
    
    # Train and predict over training and testing data
    print("\nTraining BDT...")
    print(f"Signal events: {len(signal)}")
    print(f"Background events: {len(bkg_weights)}")
    print(f"Background original N-dim reweights: {sum(bkg_weights['Ndimreweight'])}")
    print(f"dTrain events: {len(dtrain)}")
    print(f"dTest events: {len(dtest)}\n")

    # Define hyperparameter space
    space = {
        "max_depth": hp.quniform("max_depth", 3, 18, 1),
        "learning_rate": hp.uniform ("learning_rate", 0, 0.3),
        "gamma": hp.uniform ("gamma", 1, 20),
        "reg_alpha": hp.uniform("reg_alpha", 2, 20),
        "reg_lambda": hp.uniform("reg_lambda", 2, 20),
        "colsample_bytree": hp.uniform("colsample_bytree", 0.5, 1),
        "colsample_bylevel": hp.uniform("colsample_bylevel", 0.5, 1),
        "colsample_bynode": hp.uniform("colsample_bynode", 0.5, 1),
        "min_child_weight": hp.quniform("min_child_weight", 0, 10, 1),
        "n_estimators": hp.quniform("n_estimators", 50, 200, 1),
        "subsample": hp.uniform("subsample", 0, 0.95)
    }

    def gen_model(space: Optional[dict] = None) -> float:
        """
        Model setup and training in a consistent way with hyperopt.
        """

        # Initialize wandb
        """
        run = wandb.init(
            project = "H4g Event Selection BDT",
            config = space
        )
        """

        # Setup for model
        if space is not None:
            # Must use this order of parameters in the JSON
            """
            learning_rate
            gamma
            subsample
            n_estimators
            max_depth
            reg_lambda
            reg_alpha
            colsample_bytree
            colsample_bylevel
            colsample_bynode
            min_child_weight
            """

            bdt = xgboost.XGBClassifier(
                objective = "binary:logistic",
                sampling_method = "uniform",
                learning_rate = space["learning_rate"],
                gamma = space["gamma"],
                subsample = space["subsample"],
                n_estimators = int(space["n_estimators"]),
                max_depth = int(space["max_depth"]),
                reg_lambda = space["reg_lambda"],
                reg_alpha = space["reg_alpha"],
                colsample_bytree = space["colsample_bytree"],
                colsample_bylevel = space["colsample_bylevel"],
                colsample_bynode = space["colsample_bynode"],
                min_child_weight = space["min_child_weight"],
                use_label_encoder = False, # Just to remove warning since this is deprecated
                verbosity = 0 if debug is False else 3,
                n_jobs = 1
            )

            # wtrain is not only bkg so weights are not just 1.0 when specificing no reweighting
            print(f"Sum of weights = {sum(dtrain.weight)}")
            print(f"Length of wtrain = {len(dtrain.weight)}")
        else:
            # Test with no parameters
            bdt = xgboost.XGBClassifier(
                subsample = 0.8,  # Set this to 0.8 by default since 1.0 always causes overtraining.
                n_jobs = 1
            )

        # Train model
        print("\nUsing following columns for training:")
        for i, col in enumerate(dtrain.drop(columns=["category", "mass_gggg", "weight"], axis=1).columns):
            j = i + 1
            if j % 3 == 0:
                endl = "\n"
            else:
                endl = ""
            print(f"{col:22}", end=endl)
        print("\n")

        bdt.fit(
            dtrain.drop(columns=["category", "mass_gggg", "weight"], axis=1),
            dtrain.category,
            sample_weight = dtrain["weight"],
            eval_set = [(dtrain.drop(columns=["category", "mass_gggg", "weight"], axis=1), dtrain.category), (dtest.drop(columns=["category", "mass_gggg", "weight"], axis=1), dtest.category)],
            eval_metric = ["auc", "logloss"],
            early_stopping_rounds = 25,
            #callbacks = [WandbCallback()]
        )

        # Organize results for plotting
        results = {
            "train": predict_BDT(bdt, dtrain.drop(columns=["category", "mass_gggg", "weight"], axis=1)),
            "test": predict_BDT(bdt, dtest.drop(columns=["category", "mass_gggg", "weight"], axis=1)),
        }

        train_len = len(results['train'])
        test_len = len(results['test'])
        print(f"Results Train: {train_len}")
        print(f"Results Test: {test_len}")

        # Want best value of eval metric
        predictions = results["train"]
        logloss = log_loss(dtrain["category"], predictions)
        auc = roc_auc_score(dtrain["category"], predictions)
        avg_deviation_from_truth = np.sum(np.absolute(predictions - dtrain["category"].tolist())) / len(dtrain["category"])

        # Parameterize predictions to satisfy binary classification
        error_rate = np.sum(np.heaviside(predictions - 0.5, 1) != dtrain["category"]) / len(dtrain["category"])
        f1 = f1_score(dtrain["category"], np.heaviside(predictions - 0.5, 1), sample_weight=dtrain.weight)

        # Include BDT score and weights (for plotting) in results dict
        for (key, val), df, w in zip(results.items(), [dtrain, dtest], [dtrain.weight, dtest.weight]):
            results[key] = pd.concat([pd.Series(val, name="BDT_score"), df], axis=1)
            results[key] = pd.concat([w, results[key]], axis=1)

        fig = plotting.plot_test_train_superimposed(results, year = year, xgb_name = xgb_name, suppressSave = True)

        # Removed wandb connection
        """
        # Log test metrics to wandb
        params = {}
        params["Log Loss"] = logloss
        params["AUC"] = auc
        params["Avg Truth Deviation"] = avg_deviation_from_truth
        params["Error Rate"] = error_rate
        params["F1 Score"] = f1
        params["learning_rate"] = space["learning_rate"] if space is not None else 0.3
        params["gamma"] = space["gamma"] if space is not None else 0
        params["n_estimators"] = space["n_estimators"] if space is not None else 100
        params["max_depth"] = space["max_depth"] if space is not None else 6
        params["reg_lambda"] = space["reg_lambda"] if space is not None else 1
        params["reg_alpha"] = space["reg_alpha"] if space is not None else 0
        params["colsample_bytree"] = space["colsample_bytree"] if space is not None else 1
        params["colsample_bylevel"] = space["colsample_bylevel"] if space is not None else 1
        params["colsample_bynode"] = space["colsample_bynode"] if space is not None else 1
        params["min_child_weight"] = space["min_child_weight"] if space is not None else 1
        params["subsample"] = space["subsample"] if space is not None else 0.8
        params["Nsignal_train"] = len(dtrain[dtrain.category == 1])
        params["Nbkg_train"] = len(dtrain[dtrain.category == 0])
        params["Nsignal_test"] = len(dtest[dtest.category == 1])
        params["Nbkg_test"] = len(dtest[dtest.category == 0])
        params_table = wandb.Table(columns=list(params.keys()), data=[list(params.values())])

        run.log({"Parameters": params_table})
        run.log({"ROC_Curve" : wandb.plot.roc_curve(np.array(dtrain["category"]), np.stack((predictions, 1 - predictions), axis=1))})
        run.log({"BDT Score Distribution": wandb.Image(fig)})
        
        run.finish()
        """

        return {"loss": -auc, "status": STATUS_OK, "model": bdt}

    # Randomized hyperparameter tuning
    if gen is None:
        if not train_no_params:
            best = fmin(
                fn = gen_model,
                space = space,
                algo = tpe.suggest,
                max_evals = 200,
            )

            # Run optimal model from randomized sweep
            generator = gen_model(best)
            bdt = generator["model"]

        else:
            generator = gen_model()
            bdt = generator["model"]
    
    # Generate model from input parameters
    else:
        gen_path = f"{cwd}/../../scripts/bdtIO/{gen}"
        assert os.path.exists(gen_path), f"Requested parameters JSON {gen_path} does not exist."

        with open(gen_path, "r") as f:
            params = json.load(f)
        
        generator = gen_model(params)
        bdt = generator["model"]
    
    # Save model
    if save:
        if os.path.exists(modelPath):
            os.remove(modelPath)

        bdt.save_model(modelPath)
        assert os.path.exists(modelPath), "Saved model does not exist."
        print("\nSaved model to file.")

        # Save BDT parameters
        json_file = f"{cwd}/../../scripts/bdtIO/BDT_{year}_outputs.json"
        with open(json_file, "w") as f:
            #print(bdt.get_params())
            params = bdt.get_params()
            try:
                del params["sample_weight"]
            except:
                pass
            json.dump(params, f, indent=4)
        
        assert os.path.isfile(json_file), "Saved JSON file does not exist."
        print("Saved BDT parameters to file.")

    # Organize results
    results = {
        "train": predict_BDT(bdt, dtrain.drop(columns=["category", "mass_gggg", "weight"], axis=1)),
        "test": predict_BDT(bdt, dtest.drop(columns=["category", "mass_gggg", "weight"], axis=1)),
    }

    categories = {"train": dtrain.category, "test": dtest.category}

    # Remember that masss_gggg is in results for correlation matrix
    for (key, val), df in zip(results.items(), [dtrain, dtest]):
        results[key] = pd.concat([pd.Series(val, name="BDT_score"), df], axis=1)


    # ***********************************************************************------------------------------------
    # ---- Net and per-mass-point AUC (train and test) ----
    # Rebuild the same splits as prep_samples (same weight<10 cut, same even/odd rule)
    sig_all = signal[signal.weight < 10.0]
    sig_tr, sig_te = split_dataset(sig_all)
    bkg_tr, bkg_te = split_dataset(bkg_weights)

    def mass_point(df):
        # interMass = (m_reco - m_hyp) / m_gggg  ->  m_hyp
        return ((df["LeadPs_mass"] - df["LeadPs_interMass"] * df["mass_gggg"]) / 5).round() * 5

    def auc_pair(y, score, w):
        return {"auc": float(roc_auc_score(y, score)),
                "auc_weighted": float(roc_auc_score(y, score, sample_weight=w))}

    auc_summary = {}
    roc_store = {} # added for mass specific ROC
    rows = []
    for split, s, b in [("train", sig_tr, bkg_tr), ("test", sig_te, bkg_te)]:
        s_score = predict_BDT(bdt, reduce_dataframe(s))
        b_score = predict_BDT(bdt, reduce_dataframe(b))
        s_w = s["weight"].to_numpy()
        b_w = b["Ndimreweight"].to_numpy()  # global bkg rescaling cancels in a weighted AUC
        s_mpt = mass_point(s).to_numpy()

        y = np.concatenate([np.ones(len(s_score)), np.zeros(len(b_score))])
        net = auc_pair(y, np.concatenate([s_score, b_score]), np.concatenate([s_w, b_w]))
        # ---- net ROC curve (added later) ----
        fpr_n, tpr_n, _ = roc_curve(y, np.concatenate([s_score, b_score]),
                                    sample_weight=np.concatenate([s_w, b_w]))
        roc_store[(split, "net")] = (fpr_n, tpr_n)
        # -------------------------------
        
        auc_summary[split] = {"net": net, "per_mass": {}}

        for m in sorted(set(s_mpt)):
            sel = s_mpt == m
            y_m = np.concatenate([np.ones(sel.sum()), np.zeros(len(b_score))])
            sc_m = np.concatenate([s_score[sel], b_score])
            w_m = np.concatenate([s_w[sel], b_w])
            res = auc_pair(y_m, sc_m, w_m)
            fpr, tpr, _ = roc_curve(y_m, sc_m, sample_weight=w_m)
            roc_store[(split, int(m))] = (fpr, tpr)
            auc_summary[split]["per_mass"][str(int(m))] = res
            rows.append({"split": split, "mass": int(m), "n_sig": int(sel.sum()), **res})
        rows.append({"split": split, "mass": "net", "n_sig": len(s_score), **net})

    auc_dir = f"{cwd}/../../scripts/plots/BDT/{year}/{xgb_name}"
    os.makedirs(auc_dir, exist_ok=True)
    with open(f"{auc_dir}/AUC_{year}.json", "w") as f:
        json.dump(auc_summary, f, indent=4)
    auc_table = pd.DataFrame(rows)
    auc_table.to_csv(f"{auc_dir}/AUC_{year}.csv", index=False)
    print("\nAUC summary (signal per mass point vs all bkg):")
    print(auc_table.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    #-----------******************^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    np.savez(f"{auc_dir}/ROC_curves_{year}.npz",
             **{f"{sp}_{m}_{k}": v for (sp, m), (f_, t_) in roc_store.items()
                for k, v in (("fpr", f_), ("tpr", t_))})

    for split in ["train", "test"]:
        fig, ax = plt.subplots(figsize=(7, 6))
        for (sp, m), (f_, t_) in sorted(roc_store.items(), key=lambda kv: str(kv[0][1])):
            if sp != split:
                continue
            a = auc_summary[split]["net"]["auc_weighted"] if m == "net" else auc_summary[split]["per_mass"][str(m)]["auc_weighted"]
            lab = f"net (AUC={a:.3f})" if m == "net" else f"{m} GeV (AUC={a:.3f})"
            ax.plot(f_, t_, lw=2.5 if m == "net" else 1.2, ls="--" if m == "net" else "-", label=lab)
        ax.plot([0, 1], [0, 1], "k:", lw=0.8)
        ax.set_xlabel("Background efficiency"); ax.set_ylabel("Signal efficiency")
        ax.set_title(f"ROC per signal mass point ({split})"); ax.legend(fontsize=7, ncol=2)
        fig.savefig(f"{auc_dir}/ROC_per_mass_{split}_{year}.pdf"); plt.close(fig)

    # overtraining check per mass point
    for m in sorted(auc_summary["train"]["per_mass"]):
        tr = auc_summary["train"]["per_mass"][m]["auc_weighted"]
        te = auc_summary["test"]["per_mass"][m]["auc_weighted"]
        print(f"[DBG-ROC] {m} GeV: train AUC {tr:.4f} | test AUC {te:.4f} | gap {tr - te:+.4f}")
    # *****************************************************************------------------------------------------

    # Plot predictions
    plotting.plot_test_train(results, year=year, xgb_name=xgb_name)
    plotting.plot_test_train_superimposed(results, year, xgb_name=xgb_name, suppressSave = False)
    # Plot ROC curve
    plotting.plot_roc(results, f"{cwd}/../../scripts/plots/BDT/{year}/{xgb_name}/ROC_{year}.pdf", sig=len(sig_train), bkg=len(bkg_train))
    # Plot correlation matrix
    plotting.plot_corr_matrix(results, f"{cwd}/../../scripts/plots/BDT/{year}/{xgb_name}/corr_matrix_{year}.pdf")

    # Report processing time
    elapsed = time.time() - start_time
    print(f"\nElapsed training/testing time: {elapsed}s\n")

    return bdt

def predict_BDT(
    bdt: xgboost.XGBClassifier,
    dataset: pd.DataFrame
) -> pd.DataFrame:
    """
    Predict using BDT over dataset. Returns probability that input is signal rather than the inverted, native predict_proba value.
    """

    # Predictions are given as probability that input is of class "signal"
    return bdt.predict_proba(dataset)[:,1]

