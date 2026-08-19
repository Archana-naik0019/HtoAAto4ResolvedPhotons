def electron_mask(electrons):

    return (
        (electrons.pt > 25)
        &
        (abs(electrons.eta) < 2.5)
    )
