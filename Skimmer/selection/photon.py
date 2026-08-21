import awkward as ak
import numpy as np

# Detector geometry constants
GAP_BARREL_ETA = 1.4442
GAP_ENDCAP_ETA = 1.5660
MAX_ETA = 2.5


def h4g_preselection_mask(photons):
    """Generates a boolean event mask based on the 4-photon preselection sequence:

    1. Require >= 4 photons.
    2. Eta acceptance & EB/EE gap on first 4 photons.
    3. PixelSeed Veto (pixelSeed == False) on first 4 photons.
    4. pT > 12 GeV on first 4 photons.
    """
    # 1. At least 4 photons
    has_4pho = ak.num(photons, axis=1) >= 4

    # Slice first 4 photons (for events with < 4 photons, this produces empty/short records)
    pho4 = photons[:, :4]

    # 2. Eta acceptance
    abs_eta = np.abs(pho4.eta)
    eta_pass_each = (abs_eta < GAP_BARREL_ETA) | (
        (abs_eta > GAP_ENDCAP_ETA) & (abs_eta < MAX_ETA)
    )
    eta_pass = ak.all(eta_pass_each, axis=1)

    # 3. PixelSeed Veto
    pixel_pass = ak.all(pho4.pixelSeed == False, axis=1)

    # 4. pT > 15 GeV (for post-corrections)
    pt_pass = ak.all(pho4.pt > 12, axis=1)

    # Combined mask (requires all steps)
    return has_4pho & eta_pass & pixel_pass & pt_pass