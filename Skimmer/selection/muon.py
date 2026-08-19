def muon_mask(muons):

    return (
        (muons.pt > 20)
        &
        (abs(muons.eta) < 2.4)
    )
