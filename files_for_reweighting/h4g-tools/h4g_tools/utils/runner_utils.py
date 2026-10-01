import argparse
import sys
from typing import List, Optional


def get_parser_model() -> argparse.ArgumentParser:
    """
    Setup for command line arguments.
    """

    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser", conflict_handler="resolve")

    # Input arguments
    parser.add_argument("-sig", "--sigInput", help="Directory to import signal samples from.", type=str, default=None, required=False)
    parser.add_argument("-bkg", "--bkgInput", help="Directory to import background samples from.", type=str, default=None, required=False)
    parser.add_argument("-d", "--data", help="Directory to import data samples from.", type=str, default=None, required=False)

    # Test/plot/save flags
    parser.add_argument("-i", "--inPlot", help="Information in/out of plots.", default=False, action="store_true")
    parser.add_argument("-p", "--prelim", help="Draw Preliminary on plots.", default=False, action="store_true")
    parser.add_argument("-t", "--test", help="Bool for whether to test optimal number of PDF fits.", action="store_true")
    parser.add_argument("-s", "--save", help="Save the workspace to file with name...", action="store_true", required=False)

    return parser

def get_parser_BDT() -> argparse.ArgumentParser:
    """
    Setup for command line arguments relevant to running BDT.
    """

    # Use standard parser as a base
    parser = get_parser_model()

    # Arguments
    parser.add_argument("-s", "--save", help="Save outputs (model or predictions). Plotting is automatic.", action="store_true", required=False)
    parser.add_argument("-t", "--train", help="Optional argument to train BDT.", action="store_true", required=False)
    parser.add_argument("-m", "--model", help="Path of model to load/save file (within models directory).", type=str, default=None, required=False)
    parser.add_argument("-d", "--data", help="Directory to import data samples from.", type=str, default=None, required=False)
    parser.add_argument("--debug", help="Specify debugging verbosity for BDT training.", action="store_true", required=False)

    return parser

def get_parser_scan() -> argparse.ArgumentParser:
    """
    Setup for command line arguments relevant to scanning pseudoscalar masses.
    """

    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser")
    group = parser.add_mutually_exclusive_group(required=True)

    # Arguments
    parser.add_argument("-sig", "--sigInput", help="Directory to import signal samples from. Exclude mass specific directories from path.", type=str, default=None, required=False)
    parser.add_argument("-bkg", "--bkgInput", help="Directory to import background samples from.", type=str, default=None, required=False)
    parser.add_argument("-d", "--data", help="Directory to import data samples from.", type=str, default=None, required=False)
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    parser.add_argument("-xgb", "--xgb_name", type=str, required=True, help="XGBoost model name.")
    parser.add_argument("--debug", default=False, action="store_true", help="Turns on verbosity 3 for XGBoost.")
    parser.add_argument("-fp", "--from_params", help="Specify parameters from file for BDT model.", action="store_true", required=False)
    parser.add_argument("-tnp", "--train_no_params", help="Train BDT model without any parameters.", action="store_true", required=False)
    group.add_argument("-g", "--gen", action="store_true", help="Generate BDTs with scan.")
    group.add_argument("-r", "--run", action="store_true", help="Run BDTs over data with scan.")
    group.add_argument("-p", "--plot", action="store_true", help="Plot comparison of signal, background, and data.")

    return parser

def get_parser_plotting() -> argparse.ArgumentParser:
    """
    Setup for command line arguments relevant to scanning pseudoscalar masses.
    """

    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser")

    # Arguments
    parser.add_argument("-i", "--inputs", default=None, help="Input directory for signal samples.")
    parser.add_argument("-d", "--data", default=None, help="Input directory for data samples to examine m_hyp.")
    parser.add_argument("-bkg", "--background", default=None, help="Input directory for background samples to examine m_hyp.")
    parser.add_argument("-cb", "--cutBDT", default=None, type=float, help="Cut on BDT before plotting.")
    parser.add_argument("-noRW", dest="noRW", action="store_true", help="Suppress reweighting plots.")

    return parser

def get_parser_cats() -> argparse.ArgumentParser:
    """
    Setup for command line arguments for BDT categories script.
    """

    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser")

    # Input arguments
    parser.add_argument("-c", "--nCats", help="Specify number of categories to process.", default=1, type=int, required=False)
    parser.add_argument("-sm", "--smoothness", help="Sets bass (smoothness) level of SmoothSuper.", default=7.5, type=float, required=False)
    parser.add_argument("-sp", "--span", help="Sets span of SmoothSuper. Must be between 0 and 1.", default=0.02, type=float, required=False)
    parser.add_argument("--noBkgMin", help="Removes requirement of min 8 bkg sideband events.", action="store_true", required=False)
    #parser.add_argument("--use0p6", help="Smooths signal and calculates AMS based on BDT score >= 0.6", action="store_true", required=False)

    return parser

def get_parser_cat_cuts() -> argparse.ArgumentParser:
    """
    Setup for command line arguments for applying BDT categories script.
    """

    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser")

    # Input arguments
    parser.add_argument("-cats", "--categories", help="Specify number of categories to check.", default=1, type=int, required=False)

    return parser

def handle_errors(
    *args: List
) -> None:
    """
    Handles errors from argument parser in a clean way.
    """

    # Create dummy parser
    parser = argparse.ArgumentParser()
    
    errors = False
    for error in args:
        err = error[0]
        bad_cond = True

        # Combine bad conditions
        for e in error[1:]:
            bad_cond = bad_cond and e

        # Check bad conditions
        if bad_cond:
            try:
                parser.error(err)
            except:
                errors = True
    
    assert not errors, "Argument errors are present!"
