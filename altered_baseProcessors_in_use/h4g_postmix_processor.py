"""
This script is a modified version of the base.py, with only those steps that are to be implemented after the photon-preselections are applied.
"""

from __future__ import annotations

import argparse
import logging
import os
import pathlib
import re
import shutil
import warnings
from typing import Any, Dict, Optional

import awkward
import gc
import numpy
import pandas
import pyarrow
import pyarrow.parquet
import vector
import xgboost

from higgs_dna.tools.xgb_loader import load_bdt
from higgs_dna.utils.dumping_utils import diphoton_list_to_pandas, ps_list_to_pandas

logger = logging.getLogger(__name__)
vector.register_awkward()

#from coffea.nanoevents.methods import vector as _coffea_vector_methods

#awkward.behavior.update(_coffea_vector_methods.behavior)

import coffea.nanoevents.methods.candidate as _coffea_candidate

awkward.behavior.update(_coffea_candidate.behavior)

#this tah is required only while assigning the m_hyp to events in signal samples
_MASS_TAG_RE = re.compile(r"(\d+GeV)")


class H4gPostMixProcessor:
    """
    Runs the h4g-specific selection (produce_and_select_ps), BDT scoring,
    and final dump on an already-preselected / already-mixed parquet file.
    """

    def __init__(
        self,
        h4g_bdt_dir: str = "/eos/home-a/arnaik/HiggsDNA_LCG/higgs-dna-2024/higgs_dna/metaconditions/h4g/",
        bdt_year: int = 2022,
        output_location: Optional[str] = None,
        output_format: str = "parquet",
        signal: bool = False,
        mass_tag: Optional[str] = None,
    ) -> None:
        self.output_location = output_location
        self.output_format = output_format
        self.signal = signal
        self.data_kind = "data"  # mixed events carry no genWeight (treat like data)

        self.mass_tag = mass_tag

        # same prefix map HggBaseProcessor uses when flattening diphotons -> pandas
        self.prefixes = {
            "pho_lead": "lead",
            "pho_sublead": "sublead",
            "pho_subsublead": "subsublead",
            "pho_subsubsublead": "subsubsublead",
        }

        # h4g-specific object selection cuts -- copied from HggBaseProcessor.__init__
        self.lead_pt = 30.0
        self.sublead_pt = 18.0
        self.min_pt = 15.0
        self.eta_gap_low = 1.4442
        self.eta_gap_high = 1.566
        self.eta_max = 2.5
        self.mass_range_low = 110.0
        self.mass_range_high = 180.0

        #load the h4g event-selection BDT (copied from base.py) ---
        self.h4g_bdt = None
        try:
            bdt = xgboost.XGBClassifier()  # noqa: F841 (kept for parity with base.py)
            self.h4g_bdt = load_bdt(
                os.path.join(h4g_bdt_dir, f"ES_BDT_{bdt_year}.ubj"),
                classifier=True,
            )
            logger.info(f"Loaded h4g BDT: {h4g_bdt_dir}ES_BDT_{bdt_year}.ubj")
        except Exception as e:
            warnings.warn(f"Could not instantiate h4g BDT: {e}")
            self.h4g_bdt = None

    # ------------------------------------------------------------------
    # produce_and_select_ps -- copied from HggBaseProcessor
    # ------------------------------------------------------------------
    def produce_and_select_ps(
        self,
        meta: str,
        photons: awkward.Array,
        diphotons: awkward.Array,
        events: awkward.Array,
        signal: bool,
        variation: str,
        Nevents: dict,
    ) -> Any:
        #photons = photons[awkward.argsort(photons.pt, ascending=False, axis=1)] #disabled temporarily
        photons = photons[awkward.num(photons, axis=1) >= 4]

        photons["mass"] = numpy.zeros_like(photons.pt)
        photons["charge"] = numpy.zeros_like(photons.pt)

        photons = photons[awkward.num(photons, axis=1) >= 4]
        pseudos = awkward.combinations(photons[:, :4], 4, fields=["pho1", "pho2", "pho3", "pho4"])
        pseudos["pho1"] = awkward.with_name(pseudos.pho1, "PtEtaPhiMCandidate")
        pseudos["pho2"] = awkward.with_name(pseudos.pho2, "PtEtaPhiMCandidate")
        pseudos["pho3"] = awkward.with_name(pseudos.pho3, "PtEtaPhiMCandidate")
        pseudos["pho4"] = awkward.with_name(pseudos.pho4, "PtEtaPhiMCandidate")
        pseudos["mass_gggg"] = (pseudos.pho1 + pseudos.pho2 + pseudos.pho3 + pseudos.pho4).mass
        pseudos["sumPt_gggg"] = pseudos.pho1.pt + pseudos.pho2.pt + pseudos.pho3.pt + pseudos.pho4.pt

        lead_pt = pseudos.pho1.pt > self.lead_pt
        sublead_pt = pseudos.pho2.pt > self.sublead_pt
        min_pt = (pseudos.pho3.pt > self.min_pt) & (pseudos.pho4.pt > self.min_pt)
        pseudos = pseudos[awkward.firsts(lead_pt & sublead_pt & min_pt)]
        diphotons = diphotons[awkward.firsts(lead_pt & sublead_pt & min_pt)]
        events = events[awkward.firsts(lead_pt & sublead_pt & min_pt)]
        assert len(pseudos) == len(diphotons) == len(events), (len(pseudos), len(diphotons), len(events))
        Nevents.update({"pt_cuts_pseudos": len(pseudos)})

        eta_gap = (
            ((abs(pseudos.pho1.eta) < self.eta_gap_low) | (abs(pseudos.pho1.eta) > self.eta_gap_high))
            & ((abs(pseudos.pho2.eta) < self.eta_gap_low) | (abs(pseudos.pho2.eta) > self.eta_gap_high))
            & ((abs(pseudos.pho3.eta) < self.eta_gap_low) | (abs(pseudos.pho3.eta) > self.eta_gap_high))
            & ((abs(pseudos.pho4.eta) < self.eta_gap_low) | (abs(pseudos.pho4.eta) > self.eta_gap_high))
        )
        eta_max = (
            (abs(pseudos.pho1.eta) < self.eta_max)
            & (abs(pseudos.pho2.eta) < self.eta_max)
            & (abs(pseudos.pho3.eta) < self.eta_max)
            & (abs(pseudos.pho4.eta) < self.eta_max)
        )
        eta_cuts = eta_gap & eta_max
        pseudos = pseudos[awkward.firsts(eta_cuts)]
        diphotons = diphotons[awkward.firsts(eta_cuts)]
        events = events[awkward.firsts(eta_cuts)]
        assert len(pseudos) == len(diphotons) == len(events), (len(pseudos), len(diphotons), len(events))
        Nevents.update({"eta_pseudos": len(pseudos)})

        e_veto = (
            (pseudos.pho1.pixelSeed == False)
            & (pseudos.pho2.pixelSeed == False)
            & (pseudos.pho3.pixelSeed == False)
            & (pseudos.pho4.pixelSeed == False)
        )
        pseudos = pseudos[awkward.firsts(e_veto)]
        diphotons = diphotons[awkward.firsts(e_veto)]
        events = events[awkward.firsts(e_veto)]
        assert len(pseudos) == len(diphotons) == len(events), (len(pseudos), len(diphotons), len(events))
        Nevents.update({"pixelSeed_pseudos": len(pseudos)})

        mass_low = pseudos.mass_gggg > self.mass_range_low
        mass_high = pseudos.mass_gggg < self.mass_range_high
        # NOTE: mass cut left disabled temporarily, to be enabled while processing full 2024 data
        mass_cuts = awkward.ones_like(awkward.any(mass_low & mass_high, axis=1), dtype=bool)
        pseudos = pseudos[mass_cuts]
        Nevents.update({"mass": len(pseudos)})

        diphotons = diphotons[mass_cuts]
        events = events[mass_cuts]
        diphotons = diphotons[~awkward.is_none(pseudos)]
        events = events[~awkward.is_none(pseudos)]
        assert len(diphotons) == len(pseudos) == len(events), f"{len(diphotons)} {len(photons)} {len(events)}"

        assert len(pseudos) == len(pseudos[~awkward.is_none(pseudos)])

        len_ps = len(pseudos)

        
#----------------------------------------------need to alter this function slighty for pT unordered case-----------------------------
        '''
        def makePsPermutation(ps_array, perm):
            phoA = ps_array[f"pho{perm[0]}"]
            phoB = ps_array[f"pho{perm[1]}"]
            #assert awkward.all(awkward.all((phoA.pt >= phoB.pt), axis=1)), "PhoA pT must be bigger than PhoB pT!" #temporarily turning it off

            ps_perm = awkward.zip({"ps": phoA + phoB, "leading_pho": phoA, "subleading_pho": phoB})

            for field in ["ps", "leading_pho", "subleading_pho"]:
                ps_perm[field] = ps_perm[field]
                ps_perm[(field, "pt")] = ps_perm[field].pt
                ps_perm[(field, "eta")] = ps_perm[field].eta
                ps_perm[(field, "phi")] = ps_perm[field].phi
                ps_perm[(field, "mass")] = ps_perm[field].mass

            ps_perm["ps"] = awkward.with_name(ps_perm.ps, "PtEtaPhiMCandidate")
            ps_perm["leading_pho"] = awkward.with_name(ps_perm.leading_pho, "PtEtaPhiMCandidate")
            ps_perm["subleading_pho"] = awkward.with_name(ps_perm.subleading_pho, "PtEtaPhiMCandidate")

            return ps_perm
            '''
#------------------------------------------------------------------------------------------------------------------------------------
        def makePsPermutation(ps_array, perm):
            pA = awkward.firsts(ps_array[f"pho{perm[0]}"])
            pB = awkward.firsts(ps_array[f"pho{perm[1]}"])
        
            phoA_is_lead = awkward.fill_none(pA.pt >= pB.pt, False)
        
            phoA_flat = awkward.where(phoA_is_lead, pA, pB)
            phoB_flat = awkward.where(phoA_is_lead, pB, pA)
        
            ones = numpy.ones(len(phoA_flat), dtype=numpy.int64)
            phoA = awkward.with_name(awkward.unflatten(phoA_flat, ones), "PtEtaPhiMCandidate")
            phoB = awkward.with_name(awkward.unflatten(phoB_flat, ones), "PtEtaPhiMCandidate")
        
            ps_perm = awkward.zip({"ps": phoA + phoB, "leading_pho": phoA, "subleading_pho": phoB})
        
            for field in ["ps", "leading_pho", "subleading_pho"]:
                ps_perm[field] = ps_perm[field]
                ps_perm[(field, "pt")] = ps_perm[field].pt
                ps_perm[(field, "eta")] = ps_perm[field].eta
                ps_perm[(field, "phi")] = ps_perm[field].phi
                ps_perm[(field, "mass")] = ps_perm[field].mass
        
            ps_perm["ps"] = awkward.with_name(ps_perm.ps, "PtEtaPhiMCandidate")
            ps_perm["leading_pho"] = awkward.with_name(ps_perm.leading_pho, "PtEtaPhiMCandidate")
            ps_perm["subleading_pho"] = awkward.with_name(ps_perm.subleading_pho, "PtEtaPhiMCandidate")
        
            return ps_perm
#---------------------------------------------------------------------------------------------------------------------------------------

        Ps11 = makePsPermutation(pseudos, (1, 2))
        print(Ps11.ps.mass.type)   # debug line
        
        Ps12 = makePsPermutation(pseudos, (3, 4))
        dM1 = abs(Ps11.ps.mass - Ps12.ps.mass)

        Ps21 = makePsPermutation(pseudos, (1, 3))
        Ps22 = makePsPermutation(pseudos, (2, 4))
        dM2 = abs(Ps21.ps.mass - Ps22.ps.mass)

        Ps31 = makePsPermutation(pseudos, (1, 4))
        Ps32 = makePsPermutation(pseudos, (2, 3))
        dM3 = abs(Ps31.ps.mass - Ps32.ps.mass)

        Ps11[("ps", "gg_inEB")] = (Ps11.leading_pho.eta < self.eta_gap_low) & (Ps11.subleading_pho.eta < self.eta_gap_low)
        Ps12[("ps", "gg_inEB")] = (Ps12.leading_pho.eta < self.eta_gap_low) & (Ps12.subleading_pho.eta < self.eta_gap_low)
        Ps21[("ps", "gg_inEB")] = (Ps21.leading_pho.eta < self.eta_gap_low) & (Ps21.subleading_pho.eta < self.eta_gap_low)
        Ps22[("ps", "gg_inEB")] = (Ps22.leading_pho.eta < self.eta_gap_low) & (Ps22.subleading_pho.eta < self.eta_gap_low)
        Ps31[("ps", "gg_inEB")] = (Ps31.leading_pho.eta < self.eta_gap_low) & (Ps31.subleading_pho.eta < self.eta_gap_low)
        Ps32[("ps", "gg_inEB")] = (Ps32.leading_pho.eta < self.eta_gap_low) & (Ps32.subleading_pho.eta < self.eta_gap_low)

        mask11 = Ps11.ps.pt >= Ps12.ps.pt
        mask12 = Ps11.ps.pt < Ps12.ps.pt
        mask21 = Ps21.ps.pt >= Ps22.ps.pt
        mask22 = Ps21.ps.pt < Ps22.ps.pt
        mask31 = Ps31.ps.pt >= Ps32.ps.pt
        mask32 = Ps31.ps.pt < Ps32.ps.pt

       
        '''
        dM1_safe = awkward.where(numpy.isnan(dM1), numpy.float32(numpy.inf), dM1)
        dM2_safe = awkward.where(numpy.isnan(dM2), numpy.float32(numpy.inf), dM2)
        dM3_safe = awkward.where(numpy.isnan(dM3), numpy.float32(numpy.inf), dM3)

        dM_min = awkward.min(awkward.concatenate((dM1_safe, dM2_safe, dM3_safe), axis=1), axis=1)
        dM_mask1 = dM_min == dM1_safe
        dM_mask2 = (dM_min == dM2_safe) & ~dM_mask1
        dM_mask3 = (dM_min == dM3_safe) & ~dM_mask1 & ~dM_mask2
        '''
        #the  above section is temporarily replaced by the section below, to run a check---------------------

        dM1_safe = awkward.fill_none(awkward.firsts(dM1), numpy.inf)
        dM2_safe = awkward.fill_none(awkward.firsts(dM2), numpy.inf)
        dM3_safe = awkward.fill_none(awkward.firsts(dM3), numpy.inf)
        
        dM1_safe = awkward.nan_to_num(dM1_safe, nan=numpy.inf)
        dM2_safe = awkward.nan_to_num(dM2_safe, nan=numpy.inf)
        dM3_safe = awkward.nan_to_num(dM3_safe, nan=numpy.inf)
        
        # 1d arrays of shape (N,)
        dM_min = numpy.minimum(numpy.minimum(dM1_safe, dM2_safe), dM3_safe)
        
        dM_mask1 = dM_min == dM1_safe
        dM_mask2 = (dM_min == dM2_safe) & ~dM_mask1
        dM_mask3 = (dM_min == dM3_safe) & ~dM_mask1 & ~dM_mask2

        print("dM1 type:      ", dM1.type) #debug
        print("dM_mask1 type: ", dM_mask1.type) #debug
        print("masked type:   ", dM1[dM_mask1].type) #debug
        
        ones = numpy.ones(len(dM_min), dtype=numpy.int64)
        # wraping them back to 2D single-element lists using ones
        dM_mask1 = awkward.unflatten(dM_mask1, ones)
        dM_mask2 = awkward.unflatten(dM_mask2, ones)
        dM_mask3 = awkward.unflatten(dM_mask3, ones)
        #------------------------------------------------------------------------------------------

        LeadPs1 = awkward.concatenate((Ps11[mask11], Ps12[mask12]), axis=1)
        SubleadPs1 = awkward.concatenate((Ps11[mask12], Ps12[mask11]), axis=1)
        LeadPs2 = awkward.concatenate((Ps21[mask21], Ps22[mask22]), axis=1)
        SubleadPs2 = awkward.concatenate((Ps21[mask22], Ps22[mask21]), axis=1)
        LeadPs3 = awkward.concatenate((Ps31[mask31], Ps32[mask32]), axis=1)
        SubleadPs3 = awkward.concatenate((Ps31[mask32], Ps32[mask31]), axis=1)

        dM = awkward.concatenate((dM1[dM_mask1], dM2[dM_mask2], dM3[dM_mask3]), axis=1)
        LeadPs = awkward.concatenate((LeadPs1[dM_mask1], LeadPs2[dM_mask2], LeadPs3[dM_mask3]), axis=1)
        SubleadPs = awkward.concatenate((SubleadPs1[dM_mask1], SubleadPs2[dM_mask2], SubleadPs3[dM_mask3]), axis=1)

        pseudos["dM"] = dM
        pseudos["LeadPs"] = LeadPs.ps
        pseudos[("LeadPs", "leading_pho")] = LeadPs.leading_pho
        pseudos[("LeadPs", "subleading_pho")] = LeadPs.subleading_pho
        pseudos["SubleadPs"] = SubleadPs.ps
        pseudos[("SubleadPs", "leading_pho")] = SubleadPs.leading_pho
        pseudos[("SubleadPs", "subleading_pho")] = SubleadPs.subleading_pho

        assert len_ps == len(pseudos)

        for field in ["pho1", "pho2", "pho3", "pho4", "LeadPs", "SubleadPs"]:
            pseudos[field] = pseudos[field]
            pseudos[(field, "pt")] = pseudos[field].pt
            pseudos[(field, "eta")] = pseudos[field].eta
            pseudos[(field, "phi")] = pseudos[field].phi
            pseudos[(field, "mass")] = pseudos[field].mass

            if field != "LeadPs" and field != "SubleadPs":
                pseudos[(field, "pfRelIso_photon_pt")] = (
                    pseudos[field].pfRelIso03_chg
                    if hasattr(pseudos[field], "pfRelIso03_chg")
                    else pseudos[field].pfRelIso03_chg_quadratic
                ) * pseudos[field].pt

        pseudos["LeadPs"] = awkward.with_name(pseudos.LeadPs, "PtEtaPhiMCandidate")
        pseudos["SubleadPs"] = awkward.with_name(pseudos.SubleadPs, "PtEtaPhiMCandidate")
        pseudos["pT1_ma1"] = pseudos.LeadPs.leading_pho.pt / pseudos.LeadPs.mass
        pseudos["pT2_ma1"] = pseudos.LeadPs.subleading_pho.pt / pseudos.LeadPs.mass
        pseudos["pT1_ma2"] = pseudos.SubleadPs.leading_pho.pt / pseudos.SubleadPs.mass
        pseudos["pT2_ma2"] = pseudos.SubleadPs.subleading_pho.pt / pseudos.SubleadPs.mass

        pseudos["dR_aa"] = pseudos.LeadPs.delta_r(pseudos.SubleadPs)
        pseudos["dR_aa_mass_gggg"] = pseudos.dR_aa / pseudos.mass_gggg

        pseudos[("LeadPs", "dR_gg")] = pseudos.LeadPs.leading_pho.delta_r(pseudos.LeadPs.subleading_pho)
        pseudos[("SubleadPs", "dR_gg")] = pseudos.SubleadPs.leading_pho.delta_r(pseudos.SubleadPs.subleading_pho)

        pho1_Aframe = pseudos.LeadPs.leading_pho.boost(-pseudos.LeadPs.boostvec)
        pseudos["cos_ag"] = numpy.cos(pho1_Aframe.theta)

        if not signal:
            m_hyp = numpy.random.choice([float(x) for x in range(15, 65, 5)], len(pseudos.LeadPs.mass))
            assert len(m_hyp) == len(pseudos.LeadPs.mass)
        else:
            match = _MASS_TAG_RE.search(meta) if meta else None
            if match is None:
                raise ValueError(
                    f"--signal was set but no '<digits>GeV' mass tag could be resolved "
                    f"(checked --mass-tag override and the input filename; got "
                    f"meta={meta!r}). Pass --mass-tag explicitly, or name your input "
                    f"file with a 'NNGeV' pattern (e.g. 'mixed_signal_15GeV.parquet')."
                )
            mass_val = float(match.group(1).replace("GeV", ""))

            m_hyp = numpy.full_like(pseudos.LeadPs.mass, mass_val)
            assert len(m_hyp) == len(pseudos.LeadPs.mass)

        pseudos["m_hyp"] = m_hyp
        diphotons["masspoint"] = m_hyp
        pseudos["LeadPs_interMass"] = (pseudos.LeadPs.mass - m_hyp) / pseudos.mass_gggg
        pseudos["SubleadPs_interMass"] = (pseudos.SubleadPs.mass - m_hyp) / pseudos.mass_gggg
        pseudos["Ps_massDiff"] = pseudos.LeadPs.mass - pseudos.SubleadPs.mass

        return pseudos, diphotons, events, Nevents

    # ------------------------------------------------------------------
    def _predict_bdt(self, ps: awkward.Array) -> awkward.Array:
        features = pandas.DataFrame(
            {
                "pho1_mvaID": ps.pho1.mvaID.to_numpy().flatten(),
                "pho2_mvaID": ps.pho2.mvaID.to_numpy().flatten(),
                "pho3_mvaID": ps.pho3.mvaID.to_numpy().flatten(),
                "pho4_mvaID": ps.pho4.mvaID.to_numpy().flatten(),
                "LeadPs_pt": ps.LeadPs.pt.to_numpy().flatten(),
                "SubleadPs_pt": ps.SubleadPs.pt.to_numpy().flatten(),
                "dR_aa_mass_gggg": ps.dR_aa_mass_gggg.to_numpy().flatten(),
                "LeadPs_interMass": ps.LeadPs_interMass.to_numpy().flatten(),
                "SubleadPs_interMass": ps.SubleadPs_interMass.to_numpy().flatten(),
                "Ps_massDiff": ps.Ps_massDiff.to_numpy().flatten(),
                "cos_ag": ps.cos_ag.to_numpy().flatten(),
            }
        )
        results = self.h4g_bdt.predict_proba(features)
        ps["BDT_score"] = results[:, 1]
        return ps

    # ------------------------------------------------------------------
    def _resolve_mass_tag(self, input_path: str) -> Optional[str]:
        """
        Used to assign m_hyp, figures out the "NNGeV"-style tag to stamp onto the output filename.
        """
        if self.mass_tag:
            return self.mass_tag

        match = _MASS_TAG_RE.search(os.path.basename(input_path))
        if match:
            return match.group(1)

        return None

    # ------------------------------------------------------------------
    def run(
        self,
        input_path: str,
        dataset_name: str,
        variation: str = "nominal",
        chunk_rows: Optional[int] = None,
    ) -> Optional[pandas.DataFrame]:
        
        if chunk_rows is None:
            events = awkward.from_parquet(input_path)
            logger.info(f"Loaded {len(events)} mixed events from {input_path} (no chunking)")
            return self._run_chunk(events, input_path, dataset_name, variation, chunk_tag=None)

        pf = pyarrow.parquet.ParquetFile(input_path)
        total_rows = pf.metadata.num_rows
        logger.info(
            f"{input_path}: {total_rows} total rows across {pf.metadata.num_row_groups} row groups -- "
            f"processing chunked (chunk_rows={chunk_rows})"
        )

        last_df = None
        rows_seen = 0
        chunk_idx = 0
        pending_tables = []
        pending_rows = 0

        def _process_slice(table):
            nonlocal last_df, chunk_idx
            events = awkward.from_arrow(table)
            logger.info(f"[chunk {chunk_idx}] {table.num_rows} events")
            last_df = self._run_chunk(
                events, input_path, dataset_name, variation, chunk_tag=str(chunk_idx)
            )
            del events
            gc.collect()
            chunk_idx += 1

        def _drain(force=False):
            """
		This along with the above section is added to manual submission (not via run_analysis.py in HiggsDNA)
            """
            nonlocal pending_tables, pending_rows
            if not pending_tables:
                return
            table = pyarrow.concat_tables(pending_tables)
            offset = 0
            while table.num_rows - offset >= chunk_rows:
                _process_slice(table.slice(offset, chunk_rows))
                offset += chunk_rows
            remainder = table.slice(offset)
            if force and remainder.num_rows > 0:
                _process_slice(remainder)
                pending_tables, pending_rows = [], 0
            elif remainder.num_rows > 0:
                pending_tables, pending_rows = [remainder], remainder.num_rows
            else:
                pending_tables, pending_rows = [], 0
            del table

        for rg_idx in range(pf.metadata.num_row_groups):
            table = pf.read_row_group(rg_idx)
            pending_tables.append(table)
            pending_rows += table.num_rows
            rows_seen += table.num_rows
            if pending_rows >= chunk_rows:
                _drain()

        _drain(force=True)

        logger.info(f"Finished {input_path}: {rows_seen} rows across {chunk_idx} chunk(s)")
        return last_df

    # ------------------------------------------------------------------
    def _run_chunk(
        self,
        events: awkward.Array,
        input_path: str,
        dataset_name: str,
        variation: str,
        chunk_tag: Optional[str],
    ) -> Optional[pandas.DataFrame]:

        Nevents: Dict[str, int] = {"initial_events": len(events)}

        # Resolve the mass tag ONCE here, and thread the same value through to both produce_and_select_ps (m_hyp, signal mode only) and _write_output (output filename)
        resolved_mass_tag = self._resolve_mass_tag(input_path)

        photon_fields = awkward.fields(events.Photon)
        photons = awkward.zip(
            {f: events.Photon[f] for f in photon_fields},
            with_name="PtEtaPhiMCandidate",
        )
        #photons = photons[awkward.argsort(photons.pt, axis=1, ascending=False)], disabled temporarily
        has4 = awkward.num(photons, axis=1) >= 4
        photons, events = photons[has4], events[has4]
        Nevents["has_4_photons"] = len(photons)

        #rebuild `diphotons` from the flat, already chosen candidate fields ---
        diphotons = awkward.zip(
            {
                "pt": events.pt,
                "eta": events.eta,
                "phi": events.phi,
                "mass": events.mass,
                "charge": events.charge,
                "event_number": events.event_number,
                "run_number": events.run_number,
                "luminosity_block": events.luminosity_block,
                "rapidity": events.rapidity,
                "pT1_m_gg": events.pT1_m_gg,
                "pT2_m_gg": events.pT2_m_gg,
            },
            with_name="PtEtaPhiMCandidate",
        )
        assert len(diphotons) == len(photons) == len(events), f"{len(diphotons)} {len(photons)} {len(events)}"

        # h4g-specific object selections + pseudoscalar building ---
        pseudos, diphotons, events, Nevents = self.produce_and_select_ps(
            resolved_mass_tag, photons, diphotons, events, signal=self.signal, variation=variation, Nevents=Nevents
        )
        assert len(pseudos) == len(diphotons) == len(events), f"{len(pseudos)} {len(diphotons)} {len(events)}"
        Nevents["selections"] = len(pseudos)

        if len(pseudos) == 0:
            logger.info("No surviving events after h4g selections, nothing to write.")
            return None

        # h4g BDT score ---
        if self.h4g_bdt is not None:
            pseudos = self._predict_bdt(pseudos)
        else:
            warnings.warn("h4g BDT not loaded -- output will not have a BDT_score column.")

        # trimming photon collection to the 4 used, copied from HggBaseProcessor
        events = awkward.with_field(
            events,
            awkward.concatenate(
                [
                    awkward.firsts(pseudos["pho1"])[:, numpy.newaxis],
                    awkward.firsts(pseudos["pho2"])[:, numpy.newaxis],
                    awkward.firsts(pseudos["pho3"])[:, numpy.newaxis],
                    awkward.firsts(pseudos["pho4"])[:, numpy.newaxis],
                ],
                axis=1,
            ),
            "Photon",
        )

        
        diphotons["weight_central"] = awkward.ones_like(diphotons.pt)
        diphotons["weight"] = awkward.ones_like(diphotons.pt)

        diphotons = awkward.with_field(diphotons, awkward.firsts(pseudos["pho1"]), "pho_lead")
        diphotons = awkward.with_field(diphotons, awkward.firsts(pseudos["pho2"]), "pho_sublead")
        diphotons = awkward.with_field(diphotons, awkward.firsts(pseudos["pho3"]), "pho_subsublead")
        diphotons = awkward.with_field(diphotons, awkward.firsts(pseudos["pho4"]), "pho_subsubsublead")

        for k, v in Nevents.items():
            logger.info(f"{k}: {v}")

        # using existing dumping_utils ---
        df = diphoton_list_to_pandas(self, diphotons)
        df_ps = ps_list_to_pandas(pseudos)

        extra_cols = [df, df_ps]
        if "npvs" in events.fields:
            extra_cols.append(pandas.Series(events.npvs.to_list(), name="npvs"))
        if "npvsGood" in events.fields:
            extra_cols.append(pandas.Series(events.npvsGood.to_list(), name="npvsGood"))
        if "mix_cycle" in events.fields:
            extra_cols.append(pandas.Series(events.mix_cycle.to_list(), name="mix_cycle"))

        final_df = pandas.concat(extra_cols, axis=1)

        to_remove = [col for col in final_df.columns if "EELeak" in col]
        final_df.drop(to_remove, inplace=True, errors="ignore")

        del photons, diphotons, pseudos, events, df, df_ps, extra_cols

        if self.output_location is not None:
            self._write_output(final_df, input_path, dataset_name, variation, resolved_mass_tag, chunk_tag)

        return final_df

    # ------------------------------------------------------------------
    def _write_output(
        self,
        final_df: pandas.DataFrame,
        input_path: str,
        dataset_name: str,
        variation: str,
        mass_tag: Optional[str] = None,
        chunk_tag: Optional[str] = None,
    ) -> None:
        stem = pathlib.Path(input_path).stem

        if mass_tag and mass_tag not in stem:
            fname = f"{stem}_{mass_tag}_postmix.{self.output_format}"
        else:
            fname = f"{stem}_postmix.{self.output_format}"

        if chunk_tag is not None:
            base, ext = os.path.splitext(fname)
            fname = f"{base}_chunk{chunk_tag}{ext}"

        final_table = pyarrow.Table.from_pandas(final_df)

        subdirs = [dataset_name, variation]
        local_file = os.path.join(".", fname)
        destination = os.path.join(self.output_location, os.path.sep.join(subdirs), fname)

        pyarrow.parquet.write_table(final_table, local_file)
        dirname = os.path.dirname(destination)
        if not os.path.exists(dirname):
            pathlib.Path(dirname).mkdir(parents=True, exist_ok=True)
        shutil.copy(local_file, destination)
        assert os.path.isfile(destination)
        pathlib.Path(local_file).unlink()
        logger.info(f"Wrote {destination}")


def main() -> None:
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to the mixed_events_output.parquet file")
    parser.add_argument("--output-dir", required=True, help="Directory to write the final postmix parquet into")
    parser.add_argument("--dataset-name", default="h4g_mixed", help="Used for the output subdirectory")
    parser.add_argument("--variation", default="nominal")
    parser.add_argument("--bdt-dir", default="/eos/home-a/arnaik/HiggsDNA_LCG/higgs-dna-2024/higgs_dna/metaconditions/h4g/")
    parser.add_argument("--bdt-year", type=int, default=2022)
    parser.add_argument("--signal", action="store_true", help="Set if these are mixed *signal* events, not background")
    parser.add_argument(
        "--mass-tag",
        default=None,
        help=(
            "Manually set the mass-point tag (e.g. '15GeV') appended to the output "
            "filename and used for m_hyp in --signal mode. If omitted, it's "
            "auto-detected from the input filename (a '<digits>GeV' pattern); "
            "--dataset-name is not consulted for this."
        ),
    )
    parser.add_argument(
        "--chunk-rows",
        type=int,
        default=None,
        help=(
            "If set, stream the input parquet row-group by row-group instead of loading "
            "the whole file at once, grouping consecutive row groups until at least this "
            "many rows are accumulated before processing/writing that batch. Use this for "
            "combinatorially-mixed files where the true event count is far larger than "
            "expected (check with `pyarrow.parquet.ParquetFile(path).metadata.num_rows` "
            "first). Omit for the old whole-file-at-once behaviour."
        ),
    )
    args = parser.parse_args()

    processor = H4gPostMixProcessor(
        h4g_bdt_dir=args.bdt_dir,
        bdt_year=args.bdt_year,
        output_location=args.output_dir,
        signal=args.signal,
        mass_tag=args.mass_tag,
    )
    processor.run(
        args.input,
        dataset_name=args.dataset_name,
        variation=args.variation,
        chunk_rows=args.chunk_rows,
    )


if __name__ == "__main__":
    main()
