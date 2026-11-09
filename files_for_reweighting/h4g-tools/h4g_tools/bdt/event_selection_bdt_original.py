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

def split_dataset(
    dataset: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split dataset into training and testing. Use even and odd rows.
    """

    test = dataset.iloc[lambda x: x.index % 2 == 0]  # Even rows only
    train = dataset.iloc[lambda x: x.index % 2 == 1]  # Odd rows only

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

    # Prepare samples for training
    signal = pd.concat([signal, pd.Series(mass_gggg["signal"], name="mass_gggg")], axis=1)
    bkg_weights = pd.concat([bkg_weights, pd.Series(mass_gggg["background"], name="mass_gggg")], axis=1)
    sig_train, bkg_train, dtrain, dtest = prep_samples(signal, bkg_weights, year, Scales=Scales)

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

