from h4g_tools.utils.loading import loadPerMass

if __name__ == "__main__":

    print("~~~ Signal samples with trigger eff ~~~")
    sig_samples = loadPerMass("/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2024_SaS_Pileup_NFs_Material_FNUF_eVetoSF_PreselSFpostBPix_TriggerEff_updatedPresel_wBDT_19May2026/", branches=["weight"])
    sig_samples = dict(sorted(sig_samples.items()))

    print("\n~~~ Signal samples without trigger eff ~~~")
    no_trigger_samples = loadPerMass("/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2024_SaS_Pileup_NFs_Material_FNUF_eVetoSF_PreselSFpostBPix_updatedPresel_wBDT_04May2026/", branches=["weight"])
    no_trigger_samples = dict(sorted(no_trigger_samples.items()))

    """
    print("\n~~~ Signal samples with trigger bit ~~~")
    hlt_bit_samples = loadPerMass("/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2024_SaS_Pileup_NFs_Material_FNUF_eVetoSF_PreselSFpostBPix_HLTbit_updatedPresel_wBDT_19May2026", branches=["weight"])
    hlt_bit_samples = dict(sorted(hlt_bit_samples.items()))
    """

    print()
    for mass in sig_samples.keys():
        #print(mass.replace("Signal_","").replace("_"," "), "sum of weights:", f"TriggerEff = {sum(sig_samples[mass].weight): .4f}", "  No TriggerEff =", " " if "60" not in mass else "", f"{sum(no_trigger_samples[mass].weight): >6.4f}", f"HLT bit = {sum(hlt_bit_samples[mass].weight): .4f}", f"  Eff = {sum(sig_samples[mass].weight) / sum(no_trigger_samples[mass].weight) * 100: 2.2f}%", f"\t # events raw (trigger eff/no trigger eff/hlt bit): {len(sig_samples[mass])}", " " if "60" not in mass else "", f"{len(no_trigger_samples[mass])}", " " if "60" not in mass else "", f"{len(hlt_bit_samples[mass])}")

        print(mass.replace("Signal_","").replace("_"," "), "sum of weights:", f"TriggerEff = {sum(sig_samples[mass].weight): .4f}", "  No TriggerEff =", " " if "60" not in mass else "", f"{sum(no_trigger_samples[mass].weight): >6.4f}", f"  Eff = {sum(sig_samples[mass].weight) / sum(no_trigger_samples[mass].weight) * 100: 2.2f}%", f"\t # events raw (trigger eff/no trigger eff): {len(sig_samples[mass])}", " " if "60" not in mass else "", f"{len(no_trigger_samples[mass])}")
