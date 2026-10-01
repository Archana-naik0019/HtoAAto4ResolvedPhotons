import numpy as np
import pandas as pd
import ROOT


def recalc_nosmear(
    sample: pd.DataFrame,
) -> pd.DataFrame:
    """
    Recalcaulate mass (4g, a1, a2) with the unsmeared pT. There will be a bias due to selections applied to corrected photons, but this will compare exactly the same events between smearing and no smearing.
    """

    ### Make empty vectors for all objects
    # 4vectors for ma1
    pho1a1 = ROOT.Math.PtEtaPhiMVector()
    pho2a1 = ROOT.Math.PtEtaPhiMVector()

    # 4vectors for ma2
    pho1a2 = ROOT.Math.PtEtaPhiMVector()
    pho2a2 = ROOT.Math.PtEtaPhiMVector()

    # 4vectors for m4g
    pho1 = ROOT.Math.PtEtaPhiMVector()
    pho2 = ROOT.Math.PtEtaPhiMVector()
    pho3 = ROOT.Math.PtEtaPhiMVector()
    pho4 = ROOT.Math.PtEtaPhiMVector()

    ### Make pre-allocated arrays for recalculated masses (note: no copy call when pre-allocated)
    ma1 = np.zeros(len(sample))
    ma2 = np.zeros(len(sample))
    m4g = np.zeros(len(sample))
    
    ### Extract arrays for relevant variables
    # a1
    LeadPs_leading_pho_pt_raw = sample.LeadPs_leading_pho_pt_raw.to_numpy()
    LeadPs_leading_pho_eta = sample.LeadPs_leading_pho_eta.to_numpy()
    LeadPs_leading_pho_phi = sample.LeadPs_leading_pho_phi.to_numpy()
    LeadPs_leading_pho_M = np.zeros(len(sample))

    LeadPs_subleading_pho_pt_raw = sample.LeadPs_subleading_pho_pt_raw.to_numpy()
    LeadPs_subleading_pho_eta = sample.LeadPs_subleading_pho_eta.to_numpy()
    LeadPs_subleading_pho_phi = sample.LeadPs_subleading_pho_phi.to_numpy()
    LeadPs_subleading_pho_M = np.zeros(len(sample))

    # a2
    SubleadPs_leading_pho_pt_raw = sample.SubleadPs_leading_pho_pt_raw.to_numpy()
    SubleadPs_leading_pho_eta = sample.SubleadPs_leading_pho_eta.to_numpy()
    SubleadPs_leading_pho_phi = sample.SubleadPs_leading_pho_phi.to_numpy()
    SubleadPs_leading_pho_M = np.zeros(len(sample))

    SubleadPs_subleading_pho_pt_raw = sample.SubleadPs_subleading_pho_pt_raw.to_numpy()
    SubleadPs_subleading_pho_eta = sample.SubleadPs_subleading_pho_eta.to_numpy()
    SubleadPs_subleading_pho_phi = sample.SubleadPs_subleading_pho_phi.to_numpy()
    SubleadPs_subleading_pho_M = np.zeros(len(sample))

    # 4g
    pho1_pt_raw = sample.pho1_pt_raw.to_numpy()
    pho1_eta = sample.pho1_eta.to_numpy()
    pho1_phi = sample.pho1_phi
    pho1_M = np.zeros(len(sample))
    
    pho2_pt_raw = sample.pho2_pt_raw.to_numpy()
    pho2_eta = sample.pho2_eta.to_numpy()
    pho2_phi = sample.pho2_phi
    pho2_M = np.zeros(len(sample))

    pho3_pt_raw = sample.pho3_pt_raw.to_numpy()
    pho3_eta = sample.pho3_eta.to_numpy()
    pho3_phi = sample.pho3_phi
    pho3_M = np.zeros(len(sample))

    pho4_pt_raw = sample.pho4_pt_raw.to_numpy()
    pho4_eta = sample.pho4_eta.to_numpy()
    pho4_phi = sample.pho4_phi
    pho4_M = np.zeros(len(sample))

    ### Process new values for entire sample
    for idx in range(len(sample)):
        # Build "unsmeared" a1 components
        pho1a1.SetCoordinates(
            LeadPs_leading_pho_pt_raw[idx],
            LeadPs_leading_pho_eta[idx],
            LeadPs_leading_pho_phi[idx],
            LeadPs_leading_pho_M[idx],
        )
        pho2a1.SetCoordinates(
            LeadPs_subleading_pho_pt_raw[idx],
            LeadPs_subleading_pho_eta[idx],
            LeadPs_subleading_pho_phi[idx],
            LeadPs_subleading_pho_M[idx],
        )

        # Build "unsmeared" a2 components
        pho1a2.SetCoordinates(
            SubleadPs_leading_pho_pt_raw[idx],
            SubleadPs_leading_pho_eta[idx],
            SubleadPs_leading_pho_phi[idx],
            SubleadPs_leading_pho_M[idx],
        )
        pho2a2.SetCoordinates(
            SubleadPs_subleading_pho_pt_raw[idx],
            SubleadPs_subleading_pho_eta[idx],
            SubleadPs_subleading_pho_phi[idx],
            SubleadPs_subleading_pho_M[idx],
        )
        
        # Build "unsmeared" 4g components
        pho1.SetCoordinates(
            pho1_pt_raw[idx],
            pho1_eta[idx],
            pho1_phi[idx],
            pho1_M[idx],
        )
        pho2.SetCoordinates(
            pho2_pt_raw[idx],
            pho2_eta[idx],
            pho2_phi[idx],
            pho2_M[idx],
        )
        pho3.SetCoordinates(
            pho3_pt_raw[idx],
            pho3_eta[idx],
            pho3_phi[idx],
            pho3_M[idx],
        )
        pho4.SetCoordinates(
            pho4_pt_raw[idx],
            pho4_eta[idx],
            pho4_phi[idx],
            pho4_M[idx],
        )

        # Construct objects
        a1 = pho1a1 + pho2a1
        a2 = pho1a2 + pho2a2
        fourPho = pho1 + pho2 + pho3 + pho4

        # Set new values to masses
        ma1[idx] = a1.mass()
        ma2[idx] = a2.mass()
        m4g[idx] = fourPho.mass()

    # Reset masses with no smearing applied in original DataFrame (read: biased but unsmeared)
    assert len(sample) == len(ma1) == len(ma2) == len(m4g), (len(sample), len(ma1), len(ma2), len(m4g))
    new_sample = sample.copy()
    new_sample["LeadPs_mass"] = ma1
    new_sample["SubleadPs_mass"] = ma2
    new_sample["mass_gggg"] = m4g

    return new_sample
