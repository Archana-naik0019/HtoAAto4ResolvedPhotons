import numpy as np
from coffea.lumi_tools import LumiMask


def get_lumi_mask(run, luminosityBlock, json_path):
    """Applies a Golden-JSON luminosity mask for Data events."""
    if json_path is None:
        return np.ones(len(run), dtype=bool)

    lumi_mask = LumiMask(json_path)
    return lumi_mask(run, luminosityBlock)