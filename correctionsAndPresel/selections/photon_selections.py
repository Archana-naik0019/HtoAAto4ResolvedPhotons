import awkward as ak
import numpy


# photon preselection for Run3 -> take as input nAOD Photon collection and return the Photons that pass
# cuts (pt, eta, sieie, mvaID, iso... etc)
#
def photon_preselection(
    self,
    photons: ak.Array,
    events: ak.Array,
    electron_veto=True,
    revert_electron_veto=False,
    year="2023",
    IsFlag=False
) -> ak.Array:
    """
    Apply preselection cuts to photons.
    Note that these selections are applied on each photon, it is not based on the diphoton pair.
    """
    # hlt-mimicking cuts
    rho = events.Rho.fixedGridRhoAll * ak.ones_like(photons.pt)
    photon_abs_eta = numpy.abs(photons.eta)
    if year in ["2016", "2016PreVFP", "2016PostVFP", "2017", "2018"]:
        # Run 2, use standard photon preselection
        pass_phoIso_rho_corr_EB = (
            (photon_abs_eta < self.eta_rho_corr)
            & (
                photons.pfPhoIso03 - rho * self.low_eta_rho_corr
                < self.max_pho_iso_EB_low_r9
            )
        ) | (
            # should almost never happen because of the requirement of (photons.isScEtaEB) earlier, thus might be slightly redundant
            (photon_abs_eta > self.eta_rho_corr)
            & (
                photons.pfPhoIso03 - rho * self.high_eta_rho_corr
                < self.max_pho_iso_EB_low_r9
            )
        )

        pass_phoIso_rho_corr_EE = (
            (photon_abs_eta < self.eta_rho_corr)
            & (
                photons.pfPhoIso03 - rho * self.low_eta_rho_corr
                < self.max_pho_iso_EE_low_r9
            )
        ) | (
            (photon_abs_eta > self.eta_rho_corr)
            & (
                photons.pfPhoIso03 - rho * self.high_eta_rho_corr
                < self.max_pho_iso_EE_low_r9
            )
        )
    else:
        # quadratic EA corrections in Run3 : https://indico.cern.ch/event/1204277/contributions/5064356/attachments/2538496/4369369/CutBasedPhotonID_20221031.pdf
        pass_phoIso_rho_corr_EB = (
            ((photon_abs_eta > 0.0) & (photon_abs_eta < 1.0))
            & (
                photons.pfPhoIso03 - (rho * self.EA1_EB1) - (rho * rho * self.EA2_EB1)
                < self.max_pho_iso_EB_low_r9
            )
        ) | (
            ((photon_abs_eta > 1.0) & (photon_abs_eta < 1.4442))
            & (
                photons.pfPhoIso03 - (rho * self.EA1_EB2) - (rho * rho * self.EA2_EB2)
                < self.max_pho_iso_EB_low_r9
            )
        )

        pass_phoIso_rho_corr_EE = (
            (
                ((photon_abs_eta > 1.566) & (photon_abs_eta < 2.0))
                & (
                    photons.pfPhoIso03
                    - (rho * self.EA1_EE1)
                    - (rho * rho * self.EA2_EE1)
                    < self.max_pho_iso_EE_low_r9
                )
            )
            | (
                ((photon_abs_eta > 2.0) & (photon_abs_eta < 2.2))
                & (
                    photons.pfPhoIso03
                    - (rho * self.EA1_EE2)
                    - (rho * rho * self.EA2_EE2)
                    < self.max_pho_iso_EE_low_r9
                )
            )
            | (
                ((photon_abs_eta > 2.2) & (photon_abs_eta < 2.3))
                & (
                    photons.pfPhoIso03
                    - (rho * self.EA1_EE3)
                    - (rho * rho * self.EA2_EE3)
                    < self.max_pho_iso_EE_low_r9
                )
            )
            | (
                ((photon_abs_eta > 2.3) & (photon_abs_eta < 2.4))
                & (
                    photons.pfPhoIso03
                    - (rho * self.EA1_EE4)
                    - (rho * rho * self.EA2_EE4)
                    < self.max_pho_iso_EE_low_r9
                )
            )
            | (
                ((photon_abs_eta > 2.4) & (photon_abs_eta < 2.5))
                & (
                    photons.pfPhoIso03
                    - (rho * self.EA1_EE5)
                    - (rho * rho * self.EA2_EE5)
                    < self.max_pho_iso_EE_low_r9
                )
            )
        )

    isEB_high_r9 = (photons.isScEtaEB) & (photons.r9 > self.min_full5x5_r9_EB_high_r9)
    isEE_high_r9 = (photons.isScEtaEE) & (photons.r9 > self.min_full5x5_r9_EE_high_r9)
    iso = photons.trkSumPtHollowConeDR03 if hasattr(photons, "trkSumPtHollowConeDR03") else photons.pfChargedIsoPFPV  # photons.pfChargedIsoPFPV for v11, photons.trkSumPtHollowConeDR03 v12 and above
    rel_iso = photons.pfRelIso03_chg if hasattr(photons, "pfRelIso03_chg") else photons.pfRelIso03_chg_quadratic  # photons.pfRelIso03_chg for v1?, photons.pfRelIso03_chg_quadratic v12 and above
    isEB_low_r9 = (
        (photons.isScEtaEB)
        & (photons.r9 > self.min_full5x5_r9_EB_low_r9)
        & (photons.r9 < self.min_full5x5_r9_EB_high_r9)
        & (
            iso
            < self.max_trkSumPtHollowConeDR03_EB_low_r9
        )
        & (photons.sieie < self.max_sieie_EB_low_r9)
        & (pass_phoIso_rho_corr_EB)
    )
    isEE_low_r9 = (
        (photons.isScEtaEE)
        & (photons.r9 > self.min_full5x5_r9_EE_low_r9)
        & (photons.r9 < self.min_full5x5_r9_EE_high_r9)
        & (
            iso
            < self.max_trkSumPtHollowConeDR03_EE_low_r9
        )
        & (photons.sieie < self.max_sieie_EE_low_r9)
        & (pass_phoIso_rho_corr_EE)
    )

    if electron_veto:
        e_veto_cut = (photons.electronVeto == 1)
    elif revert_electron_veto:
        e_veto_cut = (photons.electronVeto == 0)
    else:
        e_veto_cut = ak.ones_like(photons.electronVeto, dtype=bool)

    if IsFlag:
        photons["PassPresel"] = (
            e_veto_cut
            & (photons.pt > self.min_pt_photon)
            & (photons.isScEtaEB | photons.isScEtaEE)
            & (photons.mvaID > self.min_mvaid)
            & (photons.hoe < self.max_hovere)
            & (
                (photons.r9 > self.min_full5x5_r9)
                | (
                    rel_iso * photons.pt < self.max_chad_iso
                )
                | (rel_iso < self.max_chad_rel_iso)
            )
            & (isEB_high_r9 | isEB_low_r9 | isEE_high_r9 | isEE_low_r9)
        )
        return photons
    else:
        return photons[
            e_veto_cut
            & (photons.pt > self.min_pt_photon)
            & (photons.isScEtaEB | photons.isScEtaEE)
            & (photons.mvaID > self.min_mvaid)
            & (photons.hoe < self.max_hovere)
            & (
                (photons.r9 > self.min_full5x5_r9)
                | (
                    rel_iso * photons.pt < self.max_chad_iso
                )
                | (rel_iso < self.max_chad_rel_iso)
            )
            & (isEB_high_r9 | isEB_low_r9 | isEE_high_r9 | isEE_low_r9)
        ]


def addRhoCorrections(
    self,
    photons: ak.Array,
    events: ak.Array,
) -> ak.Array:
    """
    Adds a new field for pfPhoIso03 with the EA rho corrections for Run 3.
    """

    rho = events.Rho.fixedGridRhoAll * ak.ones_like(photons.pt) if hasattr(events, "Rho") else events.fixedGridRhoAll * ak.ones_like(photons.pt)

    # Get absolute eta to keep cuts simple
    photon_abs_eta = numpy.abs(photons.eta)

    # EB regions (rhoCorr_EB_X is the region, mask_EB_X is the mask with only that region)
    rhoCorr_EB_0 = (photon_abs_eta > 0.0) & (photon_abs_eta <= 1.0)
    mask_EB_0 = ak.mask(photons.pfPhoIso03, rhoCorr_EB_0)
    pt_mask_EB_0 = ak.mask(photons.pt, rhoCorr_EB_0)
    EB_0 = mask_EB_0 - (rho * self.EA1_EB1) - (rho * rho * self.EA2_EB1)

    rhoCorr_EB_1 = (photon_abs_eta > 1.0) & (photon_abs_eta <= 1.4442)
    mask_EB_1 = ak.mask(photons.pfPhoIso03, rhoCorr_EB_1)
    pt_mask_EB_1 = ak.mask(photons.pt, rhoCorr_EB_1)
    EB_1 = mask_EB_1 - (rho * self.EA1_EB2) - (rho * rho * self.EA2_EB2)

    # EE regions (rhoCorr_EE_X is the region, mask_EE_X is the mask with only that region)
    rhoCorr_EE_0 = (photon_abs_eta > 1.566) & (photon_abs_eta <= 2.0)
    mask_EE_0 = ak.mask(photons.pfPhoIso03, rhoCorr_EE_0)
    pt_mask_EE_0 = ak.mask(photons.pt, rhoCorr_EE_0)
    EE_0 = mask_EE_0 - (rho * self.EA1_EE1) - (rho * rho * self.EA2_EE1)

    rhoCorr_EE_1 = (photon_abs_eta > 2.0) & (photon_abs_eta <= 2.2)
    mask_EE_1 = ak.mask(photons.pfPhoIso03, rhoCorr_EE_1)
    pt_mask_EE_1 = ak.mask(photons.pt, rhoCorr_EE_1)
    EE_1 = mask_EE_1 - (rho * self.EA1_EE2) - (rho * rho * self.EA2_EE2)

    rhoCorr_EE_2 = (photon_abs_eta > 2.2) & (photon_abs_eta <= 2.3)
    mask_EE_2 = ak.mask(photons.pfPhoIso03, rhoCorr_EE_2)
    pt_mask_EE_2 = ak.mask(photons.pt, rhoCorr_EE_2)
    EE_2 = mask_EE_2 - (rho * self.EA1_EE3) - (rho * rho * self.EA2_EE3)

    rhoCorr_EE_3 = (photon_abs_eta > 2.3) & (photon_abs_eta <= 2.4)
    mask_EE_3 = ak.mask(photons.pfPhoIso03, rhoCorr_EE_3)
    pt_mask_EE_3 = ak.mask(photons.pt, rhoCorr_EE_3)
    EE_3 = mask_EE_3 - (rho * self.EA1_EE4) - (rho * rho * self.EA2_EE4)

    rhoCorr_EE_4 = (photon_abs_eta > 2.4) & (photon_abs_eta <= 2.5)
    mask_EE_4 = ak.mask(photons.pfPhoIso03, rhoCorr_EE_4)
    pt_mask_EE_4 = ak.mask(photons.pt, rhoCorr_EE_4)
    EE_4 = mask_EE_4 - (rho * self.EA1_EE5) - (rho * rho * self.EA2_EE5)

    EB_EE_else = ak.mask(photons.pfPhoIso03, (~rhoCorr_EB_0) & (~rhoCorr_EB_1) & (~rhoCorr_EE_0) & (~rhoCorr_EE_1) & (~rhoCorr_EE_2) & (~rhoCorr_EE_3) & (~rhoCorr_EE_4))
    pt_EB_EE_else = ak.mask(photons.pt, (~rhoCorr_EB_0) & (~rhoCorr_EB_1) & (~rhoCorr_EE_0) & (~rhoCorr_EE_1) & (~rhoCorr_EE_2) & (~rhoCorr_EE_3) & (~rhoCorr_EE_4))

    # Apply corrections to appropriate regions and remove extra Nones to get back original shape/order
    pfPhoIso03_rhoCorrected = ak.concatenate([EB_0, EB_1, EE_0, EE_1, EE_2, EE_3, EE_4, EB_EE_else], axis=1)
    pt_correction_shuffled = ak.concatenate([pt_mask_EB_0, pt_mask_EB_1, pt_mask_EE_0, pt_mask_EE_1, pt_mask_EE_2, pt_mask_EE_3, pt_mask_EE_4, pt_EB_EE_else], axis=1)
    pfPhoIso03_rhoCorrected = pfPhoIso03_rhoCorrected[~ak.is_none(pfPhoIso03_rhoCorrected, axis=1)]
    pt_correction_shuffled = pt_correction_shuffled[~ak.is_none(pt_correction_shuffled, axis=1)]

    # Sort the corrections to recover the original order. Adding above can alter the order of the corrections relative to their respective photon!
    pfPhoIso03_rhoCorrected = pfPhoIso03_rhoCorrected[ak.argsort(pt_correction_shuffled, axis=1, ascending=False)]

    assert len(photons.pfPhoIso03) == len(pfPhoIso03_rhoCorrected)
    assert ak.all(ak.num(photons.pfPhoIso03, axis=1) == ak.num(pfPhoIso03_rhoCorrected, axis=1))
    photons["pfPhoIso03_rhoCorrected"] = pfPhoIso03_rhoCorrected

    return photons


def photon_preselection_h4g(
    self, photons: ak.Array, events: ak.Array, year="2022", Nevents: dict = None
) -> ak.Array:

    Nevents.update({"2photons": len(photons)})

    # Apply pre-selections to only the 1st/2nd RECO photon
    ###   HLT-mimicking cuts   ###
    # pT cuts - photons are already sorted by pT by construction
    pt_cuts = (photons[:,0].pt > self.min_lead_pho_pt) & (photons[:,1].pt > self.min_sublead_pho_pt)

    photons = photons[pt_cuts]
    events = events[pt_cuts]
    Nevents.update({"pt_cuts": len(photons)})

    # H/E
    hoe_eb1 = (photons[:,0].hoe < self.hoe_eb) & (numpy.abs(photons[:,0].eta) < self.gap_barrel_eta)
    hoe_ee1 = (photons[:,0].hoe < self.hoe_ee) & (numpy.abs(photons[:,0].eta) > self.gap_endcap_eta) & (numpy.abs(photons[:,0].eta) < self.max_eta)
    hoe_eb2 = (photons[:,1].hoe < self.hoe_eb) & (numpy.abs(photons[:,1].eta) < self.gap_barrel_eta)
    hoe_ee2 = (photons[:,1].hoe < self.hoe_ee) & (numpy.abs(photons[:,1].eta) > self.gap_endcap_eta) & (numpy.abs(photons[:,1].eta) < self.max_eta)

    photons = photons[(hoe_eb1 | hoe_ee1) & (hoe_eb2 | hoe_ee2)]
    events = events[(hoe_eb1 | hoe_ee1) & (hoe_eb2 | hoe_ee2)]
    Nevents.update({"hoe": len(photons)})

    # R9 in eta regions
    r9_eb_region1 = (numpy.abs(photons[:,0].eta) < self.gap_barrel_eta) & (photons[:,0].r9 > self.hlt_r9_eb)
    r9_ee_region1 = (numpy.abs(photons[:,0].eta) > self.gap_endcap_eta) & (numpy.abs(photons[:,0].eta) < self.max_eta) & (photons[:,0].r9 > self.hlt_r9_ee)
    r9_eb_region2 = (numpy.abs(photons[:,1].eta) < self.gap_barrel_eta) & (photons[:,1].r9 > self.hlt_r9_eb)
    r9_ee_region2 = (numpy.abs(photons[:,1].eta) > self.gap_endcap_eta) & (numpy.abs(photons[:,1].eta) < self.max_eta) & (photons[:,1].r9 > self.hlt_r9_ee)

    photons = photons[(r9_eb_region1 | r9_ee_region1) & (r9_eb_region2 | r9_ee_region2)]
    events = events[(r9_eb_region1 | r9_ee_region1) & (r9_eb_region2 | r9_ee_region2)]
    Nevents.update({"r9": len(photons)})

    # Sigma_ieie
    sigma_ieie_eb1 = (photons[:,0].sieie < self.sigma_ieie_eb) & (numpy.abs(photons[:,0].eta) < self.gap_barrel_eta)
    sigma_ieie_ee1 = (photons[:,0].sieie < self.sigma_ieie_ee) & (numpy.abs(photons[:,0].eta) > self.gap_endcap_eta) & (numpy.abs(photons[:,0].eta) < self.max_eta)
    sigma_ieie_eb2 = (photons[:,1].sieie < self.sigma_ieie_eb) & (numpy.abs(photons[:,1].eta) < self.gap_barrel_eta)
    sigma_ieie_ee2 = (photons[:,1].sieie < self.sigma_ieie_ee) & (numpy.abs(photons[:,1].eta) > self.gap_endcap_eta) & (numpy.abs(photons[:,1].eta) < self.max_eta)

    photons = photons[(sigma_ieie_eb1 | sigma_ieie_ee1) & (sigma_ieie_eb2 | sigma_ieie_ee2)]
    events = events[(sigma_ieie_eb1 | sigma_ieie_ee1) & (sigma_ieie_eb2 | sigma_ieie_ee2)]
    Nevents.update({"sigma_ieie": len(photons)})

    # pfPhoIso
    # quadratic EA corrections in Run3 : https://indico.cern.ch/event/1204277/contributions/5064356/attachments/2538496/4369369/CutBasedPhotonID_20221031.pdf
    pass_phoIso_rho_corr_EB1 = (photons[:,0].pfPhoIso03_rhoCorrected < self.pfPhoIso_eb) & (numpy.abs(photons[:,0].eta) < self.gap_barrel_eta)
    pass_phoIso_rho_corr_EE1 = (photons[:,0].pfPhoIso03_rhoCorrected < self.pfPhoIso_ee) & (numpy.abs(photons[:,0].eta) > self.gap_endcap_eta) & (numpy.abs(photons[:,0].eta) < self.max_eta)
    pass_phoIso_rho_corr_EB2 = (photons[:,1].pfPhoIso03_rhoCorrected < self.pfPhoIso_eb) & (numpy.abs(photons[:,1].eta) < self.gap_barrel_eta)
    pass_phoIso_rho_corr_EE2 = (photons[:,1].pfPhoIso03_rhoCorrected < self.pfPhoIso_ee) & (numpy.abs(photons[:,1].eta) > self.gap_endcap_eta) & (numpy.abs(photons[:,1].eta) < self.max_eta)

    photons = photons[(pass_phoIso_rho_corr_EB1 | pass_phoIso_rho_corr_EE1) & (pass_phoIso_rho_corr_EB2 | pass_phoIso_rho_corr_EE2)]
    events = events[(pass_phoIso_rho_corr_EB1 | pass_phoIso_rho_corr_EE1) & (pass_phoIso_rho_corr_EB2 | pass_phoIso_rho_corr_EE2)]
    Nevents.update({"photonIso": len(photons)})

    # TrackerIso
    '''
    trackerIso_eb1 = ((photons[:,0].pfChargedIsoPFPV if hasattr(photons, "pfChargedIsoPFPV") else photons[:,0].trkSumPtHollowConeDR03) < self.trackerIso_eb) & (numpy.abs(photons[:,0].eta) < self.gap_barrel_eta)
    trackerIso_ee1 = ((photons[:,0].pfChargedIsoPFPV if hasattr(photons, "pfChargedIsoPFPV") else photons[:,0].trkSumPtHollowConeDR03) < self.trackerIso_ee) & (numpy.abs(photons[:,0].eta) > self.gap_endcap_eta) & (numpy.abs(photons[:,0].eta) < self.max_eta)
    trackerIso_eb2 = ((photons[:,1].pfChargedIsoPFPV if hasattr(photons, "pfChargedIsoPFPV") else photons[:,1].trkSumPtHollowConeDR03) < self.trackerIso_eb) & (numpy.abs(photons[:,1].eta) < self.gap_barrel_eta)
    trackerIso_ee2 = ((photons[:,1].pfChargedIsoPFPV if hasattr(photons, "pfChargedIsoPFPV") else photons[:,1].trkSumPtHollowConeDR03) < self.trackerIso_ee) & (numpy.abs(photons[:,1].eta) > self.gap_endcap_eta) & (numpy.abs(photons[:,1].eta) < self.max_eta)
    '''

    trackerIso_eb1 = (photons[:,0].trkSumPtHollowConeDR03 < self.trackerIso_eb) & (numpy.abs(photons[:,0].eta) < self.gap_barrel_eta)
    trackerIso_ee1 = (photons[:,0].trkSumPtHollowConeDR03 < self.trackerIso_ee) & (numpy.abs(photons[:,0].eta) > self.gap_endcap_eta) & (numpy.abs(photons[:,0].eta) < self.max_eta)
    trackerIso_eb2 = (photons[:,1].trkSumPtHollowConeDR03 < self.trackerIso_eb) & (numpy.abs(photons[:,1].eta) < self.gap_barrel_eta)
    trackerIso_ee2 = (photons[:,1].trkSumPtHollowConeDR03 < self.trackerIso_ee) & (numpy.abs(photons[:,1].eta) > self.gap_endcap_eta) & (numpy.abs(photons[:,1].eta) < self.max_eta)

    photons = photons[(trackerIso_eb1 | trackerIso_ee1) & (trackerIso_eb2 | trackerIso_ee2)]
    events = events[(trackerIso_eb1 | trackerIso_ee1) & (trackerIso_eb2 | trackerIso_ee2)]
    Nevents.update({"trackerIso": len(photons)})

    # Tracker fiducial region
    eta = ((numpy.abs(photons[:,0].eta) < self.gap_barrel_eta) | ((numpy.abs(photons[:,0].eta) > self.gap_endcap_eta) & (numpy.abs(photons[:,0].eta) < self.max_eta))) & ((numpy.abs(photons[:,1].eta) < self.gap_barrel_eta) | ((numpy.abs(photons[:,1].eta) > self.gap_endcap_eta) & (numpy.abs(photons[:,1].eta) < self.max_eta)))

    photons = photons[eta]
    events = events[eta]
    Nevents.update({"eta": len(photons)})

    # Electron veto (use pixelSeed!)
    e_veto = (photons[:,0].pixelSeed == False) & (photons[:,1].pixelSeed == False)  # should be False to reject any pixelSeed (i.e. is a photon)

    photons = photons[e_veto]
    events = events[e_veto]
    Nevents.update({"pixelSeed": len(photons)})

    # Multiple option isolation cuts - similar here as with the hlt_all cuts
    photons.pfRelIso03_chg if hasattr(photons, "pfRelIso03_chg") else photons.pfRelIso03_chg_quadratic
    iso1 = (photons[:,0].r9 > self.iso_min_full5x5_r9) | ((photons[:,0].pfRelIso03_chg if hasattr(photons, "pfRelIso03_chg") else photons[:,0].pfRelIso03_chg_quadratic) < self.iso_max_chad) | (((photons[:,0].pfRelIso03_chg if hasattr(photons, "pfRelIso03_chg") else photons[:,0].pfRelIso03_chg_quadratic) * photons[:,0].pt) < self.iso_max_chad_rel)
    iso2 = (photons[:,1].r9 > self.iso_min_full5x5_r9) | ((photons[:,1].pfRelIso03_chg if hasattr(photons, "pfRelIso03_chg") else photons[:,1].pfRelIso03_chg_quadratic) < self.iso_max_chad) | (((photons[:,1].pfRelIso03_chg if hasattr(photons, "pfRelIso03_chg") else photons[:,1].pfRelIso03_chg_quadratic) * photons[:,1].pt) < self.iso_max_chad_rel)
    iso = iso1 & iso2

    photons = photons[iso]
    events = events[iso]
    Nevents.update({"miniAOD": len(photons)})

    # This method of doing each cut individually is what Tanay and I both do! Also, it yields a similar efficiency to Run 2!
    # This matches with HLT if we were to apply them simultaneously. This logic is consistent.
    photons = photons[ak.num(photons, axis=1) >= 2]
    events = events[ak.num(photons, axis=1) >= 2]
    photons = photons[~ak.is_none(photons)]
    events = events[~ak.is_none(events)]

    return photons, events, Nevents
