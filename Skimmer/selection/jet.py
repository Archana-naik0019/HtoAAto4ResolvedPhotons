def jet_mask(jets):

    return (
        (jets.pt > 15)
        &
        (abs(jets.eta) < 2.4)
    )
